/** @brief This class represents a "layer" of analysis. It calculates the
 *  cost of the atomic units of this layer, then folds them into a problem
 * formulation for the next layer **/
#include "folding.h"
#include <memory>
#include <string>

#include <isl/aff.h>
#include <isl/map.h>
#include <isl/polynomial.h>
#include <isl/set.h>

#include <barvinok/isl.h>
#include <barvinok/barvinok.h>
#include <barvinok/polylib.h>

/// @brief Strings representing the src and dst datum holds/requests in ISL.
struct binding_struct
{
    const std::string srcs;
    const std::string dsts;
};
typedef std::unique_ptr<binding_struct> binding;
/** @brief Defines the struct that comprises the result of folding and the unique
  * ptr to it that represents what is returned by fold. */
struct fold_struct
{
    const long cost;
    const std::string folded_repr;
};
typedef std::unique_ptr<fold_struct> fold_result;

class BranchTwig
{
    public:
        /// @brief The cost formula of the folding step for this layer.
        const std::string crease_costs;
        /** @brief The folding action to a multicastable representation after
          * calculating the cost of folding. */
        const std::string fold_formula;
        /// @brief The cost formula of the multicasting step for this layer.
        const std::string multicast_costs;
        /// @brief The collapse formulation for the next layer.
        const std::string collapse_formula;
        /// @brief The context the layer is in.
        isl_ctx *const ctx;
    public:
        /** @brief Constructs a layer from a cost formulation and a folding
          * formulation.
          * 
          * @param crease_costs The cost formula of unmulticastable datum as 
          * an ISL string. Goes under the assumption the input is either the 
          * starting geometry architecture or the output of a previous layer.
          * @param fold_formula The isl_map in a string representation that
          * projects away the unmulticastable portions of the path of a datum to
          * a dst. This is what we refer to as "folding".
          * @param multicast_formula The cost formula for multicasting a datum
          * as an ISL string. Goes under the assumption that the input is of
          * the form of this Layer's ISL representation after folding.
          * @param ctx The context the layer is in.
          */
        BranchTwig(
            const std::string& crease_costs, const std::string& fold_formula,
            const std::string& multicast_costs, const std::string& collapse_formula, 
            isl_ctx* ctx
        ):
        crease_costs(crease_costs), fold_formula(fold_formula),
        multicast_costs(multicast_costs), collapse_formula(collapse_formula), ctx(ctx) {}

        /** @brief Calculates the cost of the atomic units of this layer, then
          * returns a struct of the binding cost at layer and the next layer's
          * abstraction as a string. 
          * 
          * @param s_srcs The sources of the bindings at this layer as an ISL string.
          * @param s_dsts The destinations of the bindings at this layer as an ISL string.
          * 
          * @return A struct of the binding abstraction for the next layer.
          * */
        void evaluate(const std::string& s_srcs, const std::string& s_dsts)
        {
            // Folds the destinations onto their connected trunk.
            const fold_result fold_res = this->fold(s_dsts);
            std::cout << "Crease Cost: " << fold_res->cost << std::endl;
            std::cout << "Folded: " << fold_res->folded_repr << std::endl;

            // Calculates the cost to every folded node per datum.
            const long casting_cost = this->multicast(fold_res->folded_repr);
            std::cout << "Casting Cost: " << casting_cost << std::endl;

            // Calculates the total cost of the layer.
            const long total_cost = fold_res->cost + casting_cost;

            // Calculates the requests that are not satisfied by the layer.
            ///@todo Collapse the folded destinations into the next layer.
            const binding collapsed = this->collapse(s_srcs, s_dsts);
            std::cout << "Collapsed: " << collapsed->srcs << std::endl;
            std::cout << "Missing: " << collapsed->dsts << std::endl;
        }

        /// @brief Wraps evaluate by accepting the binding as a struct.
        void inline evaluate(const binding& b)
        {
            this->evaluate(b->srcs, b->dsts);
        }
    private:
        /** @brief Folds the destinations onto their connected trunk. 
         * 
         * @param s_dsts The destinations to fold as an ISL string.
         * @return A unique_ptr to a struct holding costs of the folding step and
         * the folded destinations as an ISL string.
         */
        fold_result fold(const std::string& dsts)
        {
            // Reads the dsts into isl format.
            isl_map *p_dsts = isl_map_read_from_str(ctx, dsts.c_str());
            
            /// @note Gets the total cost of the folded dsts.
            // Returns { [id, x, y] -> number_of_data}
            isl_pw_qpolynomial *p_card = isl_map_card(isl_map_copy(p_dsts));
            // Calculates the cost per datum per dst cast from the trunk.
            isl_pw_qpolynomial *p_fold_cost = isl_pw_qpolynomial_read_from_str(ctx, this->crease_costs.c_str());
            // Calculates the cost per dst cast from the trunk.
            isl_pw_qpolynomial *p_cost_at_dst = isl_pw_qpolynomial_mul(p_card, p_fold_cost);
            // Calculates the cost to cast all data from the trunk.
            isl_pw_qpolynomial *p_total_cost = isl_pw_qpolynomial_sum(p_cost_at_dst);
            // Reads the value from p_total_cost.
            isl_val *v_total_cost = isl_pw_qpolynomial_eval(p_total_cost, isl_point_zero(isl_pw_qpolynomial_get_domain_space(p_total_cost)));
            // Initializes the variable storing the cost of folding.
            long fold_cost = isl_val_get_num_si(v_total_cost);
            // Frees the values.
            isl_val_free(v_total_cost);

            /// @note Folds p_dsts onto the trunk according to the fold formula.
            // Reads the fold formula into isl format.
            isl_map *p_fold = isl_map_read_from_str(ctx, this->fold_formula.c_str());
            // Converts dsts->data to data->dsts
            isl_map *p_data_to_dsts = isl_map_reverse(p_dsts);
            // Folds the dsts onto the trunk.
            isl_map *p_folded = isl_map_apply_range(p_data_to_dsts, p_fold);
            // Gets the largest y value per datum.
            std::string all_after = "{ [id, y] -> [id, y'] : y' > y }";
            isl_map *p_all_after = isl_map_read_from_str(ctx, all_after.c_str());
            isl_map *p_max_y = isl_map_apply_range(isl_map_copy(p_folded), p_all_after);
            isl_map *p_folded_condensed = isl_map_subtract(p_folded, p_max_y);
            // Converts p_folded_condensed to a string.
            std::string s_folded_condensed = isl_map_to_str(p_folded_condensed);
            // Frees the maps.
            isl_map_free(p_folded_condensed);

            // Initializes the struct ptr to return.
            fold_result result = fold_result(new fold_struct{fold_cost, s_folded_condensed});

            return result;
        }

        /** @brief Calculates the the cost to every folded node per datum.
         * 
         * @param folded_geometry The folded destinations from this->fold(*).
         * @return The cost of multicasting to the folded destinations.
         */
        long multicast(const std::string& s_folded_bindings)
        {
            // Reads the folded destinations into isl format.
            isl_map *p_folded_bindings = isl_map_read_from_str(ctx, s_folded_bindings.c_str());
            // Reads the multicast cost formulation into isl format.
            isl_pw_qpolynomial *p_cast_cost = isl_pw_qpolynomial_read_from_str(ctx, this->multicast_costs.c_str());

            /** @note Calculates the cost of multicasting to the folded dsts
             * according to architecture spec. */
            // Applies the cost formulation to the folded dsts.
            isl_pw_qpolynomial *p_cost_applied = isl_map_apply_pw_qpolynomial(p_folded_bindings, p_cast_cost);
            // Sums all the costs.
            isl_pw_qpolynomial *p_total_cost = isl_pw_qpolynomial_sum(p_cost_applied);
            // Evaluates the cost.
            isl_val *v_total_cost = isl_pw_qpolynomial_eval(p_total_cost, isl_point_zero(isl_pw_qpolynomial_get_domain_space(p_total_cost)));
            // Initializes the return value to the cost of multicasting.
            long cost = isl_val_get_d(v_total_cost);
            // Frees the values.
            isl_val_free(v_total_cost);

            return cost;
        }

        /**
         * @brief Identifies the requests that are not satisfied by the layer 
         * and passes them as dsts to the next layer. Also calculates the cost 
         * of getting all srcs to a position accessible by the next layer.
         * 
         * @return The collapsed binding abstraction for the next layer.
         */
        binding collapse(const std::string& s_srcs, const std::string& s_dsts)
        {
            // Reads the srcs and dsts into isl format.
            isl_map *p_srcs = isl_map_read_from_str(ctx, s_srcs.c_str());
            isl_map *p_dsts = isl_map_read_from_str(ctx, s_dsts.c_str());
            // Reads the collapse formula into isl format.
            isl_map *p_collapse = isl_map_read_from_str(ctx, this->collapse_formula.c_str());

            // Collapses all dst requests to the same format as their SRCs.
            isl_map *p_collapsed = isl_map_apply_range(p_collapse, p_dsts);
            // Calculates the requests that are not satisfied by the layer.
            isl_map *p_missing = isl_map_subtract(p_collapsed, p_srcs);
            // Converts p_missing to a string.
            std::string s_missing = isl_map_to_str(p_missing);
            // Frees the maps.
            isl_map_free(p_missing);

            /** @note Initializes the collapsed binding abstraction for the next 
              * layer. */
            binding collapsed = binding(new binding_struct{s_srcs, s_missing});
            return collapsed;
        }

};

int main(int argc, char* argv[])
{
    // Creates an isl context.
    isl_ctx *ctx = isl_ctx_alloc();
    std::string srcs = R"SRC(
        {[id] -> [data] : id = 0 and 0 <= data <= 1}
    )SRC";
    // isl_map *src = isl_map_read_from_str(ctx, srcs.c_str());
    std::string data = R"DST(
        {[id, x, y] -> [data]: id = 0 and (-1 = x or x = 1) and 0 <= y <= 1 and data = y}
    )DST";
    binding test_case = binding(new binding_struct{srcs, data});

    // isl_map* id_to_all_x_y = isl_map_read_from_str(
    //     ctx,
    //     "{ [id] -> [id, x, y] }"
    // );
    // isl_map* id_to_all_dst_data = isl_map_apply_range(id_to_all_x_y, dst);
    // isl_map* id_to_missing_data = isl_map_subtract(id_to_all_dst_data, src);

    std::string crease_costs = "{ [id, x, y] -> x: x >= 0; [id, x, y] -> -x: x < 0 }";
    std::string fold_formula = "{ [id, x, y] -> [id, y] }";
    std::string multicast_costs = "{ [id, y] -> y }";
    std::string collapse_formula = "{ [id] -> [id, x, y] }";

    BranchTwig test = BranchTwig(crease_costs, fold_formula, multicast_costs, collapse_formula, ctx);
    std::cout << "Evaluating..." << std::endl;
    test.evaluate(test_case);

    ///@note Frees ctx to check for memory leaks through ISL.
    isl_ctx_free(ctx);
}