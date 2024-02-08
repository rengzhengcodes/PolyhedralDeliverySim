/** 
 * @brief This class represents a "layer" of analysis. It calculates the
 * cost of the atomic units of this layer, then folds them into a problem
 * formulation for the next layer.
 */
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

class BranchTwig
{
    public:
        /// @brief The cost formula of the folding step for this layer.
        const std::string crease_costs;
        /** 
         * @brief The folding action to a multicastable representation after
         * calculating the cost of folding. 
         */
        const std::string fold_formula;
        /// @brief The cost formula of the multicasting step for this layer.
        const std::string multicast_costs;
        /// @brief The src collapse formulation for the next layer.
        const std::string src_collapser;
        /// @brief The dst collapse formulation for the next layer.
        const std::string dst_collapser;
        /// @brief The context the layer is in.
        isl_ctx *const ctx;
    public:
        /** 
         * @brief Constructs a layer from a cost formulation and a folding
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
         * @param collapse_formula The collapse formulation for how to translate
         * the srcs and dsts of this layer to the inputs that work with next
         * layer.
         * @param ctx The context the layer is in.
         */
        BranchTwig(
            const std::string& crease_costs, const std::string& fold_formula,
            const std::string& multicast_costs, const collapse& collapse_formulas, 
            isl_ctx *const ctx
        ):
        crease_costs(crease_costs), fold_formula(fold_formula), multicast_costs(multicast_costs), 
        src_collapser(collapse_formulas->src_collapser), dst_collapser(collapse_formulas->dst_collapser),
        ctx(ctx) {}

        /** 
         * @brief Calculates the cost of the atomic units of this layer, then
         * returns a struct of the binding cost at layer and the next layer's
         * abstraction as a string. 
         * 
         * @param s_srcs The sources of the bindings at this layer as an ISL string.
         * @param s_dsts The destinations of the bindings at this layer as an ISL string.
         * 
         * @return A struct of the binding abstraction for the next layer.
         */
        binding evaluate(const std::string& s_srcs, const std::string& s_dsts)
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
            binding collapsed = this->collapse(s_srcs, s_dsts);
            return collapsed;
        }

        /// @brief Wraps evaluate by accepting the binding as a struct.
        binding inline evaluate(const binding& b)
        {
            return this->evaluate(b->srcs, b->dsts);
        }
    private:
        /** 
         * @brief Folds the destinations onto their connected trunk. 
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
            p_folded = isl_map_reverse(p_folded);
            std::cout << "P_Folded: " << isl_map_to_str(p_folded) << std::endl;
            // Gets the largest y value per datum.
            /// @todo Functionalize this.
            std::string all_after = "{ [id, y] -> [id, y'] : y' > y }";
            isl_map *p_all_after = isl_map_read_from_str(ctx, all_after.c_str());
            isl_map *p_max_y = isl_map_apply_range(p_all_after, isl_map_copy(p_folded));
            std::cout << "Max Y: " << isl_map_to_str(p_max_y) << std::endl;
            isl_map *p_folded_condensed = isl_map_subtract(p_folded, p_max_y);
            p_folded_condensed = isl_map_reverse(p_folded_condensed);
            // Converts p_folded_condensed to a string.
            std::string s_folded_condensed = isl_map_to_str(p_folded_condensed);
            // Frees the maps.
            isl_map_free(p_folded_condensed);

            // Initializes the struct ptr to return.
            fold_result result = fold_result(new fold_struct{fold_cost, s_folded_condensed});

            return result;
        }

        /** 
         * @brief Calculates the the cost to every folded node per datum.
         * 
         * @param folded_geometry The folded destinations from this->fold(*).
         * @return The cost of multicasting to the folded destinations.
         */
        long multicast(const std::string& s_folded_bindings)
        {
            // Reads the folded destinations into isl format.
            isl_map *p_folded_bindings = isl_map_read_from_str(ctx, s_folded_bindings.c_str());
            std::cout << s_folded_bindings << std::endl;
            // Reads the multicast cost formulation into isl format.
            isl_pw_qpolynomial *p_cast_cost = isl_pw_qpolynomial_read_from_str(ctx, this->multicast_costs.c_str());

            /** 
             * @note Calculates the cost of multicasting to the folded dsts
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
         * @param s_srcs The sources of the bindings at this layer as an ISL string.
         * @param s_dsts The destinations of the bindings at this layer as an ISL string.
         * 
         * @return The collapsed binding abstraction for the next layer.
         */
        binding collapse(const std::string& s_srcs, const std::string& s_dsts)
        {
            // Reads the srcs and dsts into isl format.
            isl_map *p_srcs = isl_map_read_from_str(ctx, s_srcs.c_str());
            isl_map *p_dsts = isl_map_read_from_str(ctx, s_dsts.c_str());
            // Reads the collapse formulas into isl format.
            isl_map *p_collapse_dsts = isl_map_read_from_str(
                ctx, this->dst_collapser.c_str()
            );
            isl_map *p_collapse_srcs = isl_map_read_from_str(
                ctx, this->src_collapser.c_str()
            );

            // Collapses all src requests.
            isl_map *p_collapsed_srcs = isl_map_apply_range(p_collapse_srcs, p_srcs);
            // Converts p_collapsed_srcs to a string.
            std::string s_collapsed_srcs = isl_map_to_str(p_collapsed_srcs);
            // Collapses all dst requests to the same format as their SRCs.
            isl_map *p_collapsed_dsts = isl_map_apply_range(p_collapse_dsts, p_dsts);

            // Calculates the requests that are not satisfied by the layer.
            isl_map *p_missing_data = isl_map_subtract(p_collapsed_dsts, p_collapsed_srcs);
            // Converts p_missing to a string.
            std::string s_missing_data = isl_map_to_str(p_missing_data);
            // Frees the maps.
            isl_map_free(p_missing_data);

            // Initializes the collapsed binding abstraction for the next layer.
            binding collapsed = binding(new binding_struct{s_collapsed_srcs, s_missing_data});
            return collapsed;
        }
};

class BranchTrunk
{
    public:
        /// @brief The cost formula of the folding step for this layer.
        const std::string crease_costs;
        /** 
         * @brief The folding action to a multicastable representation after
         * calculating the cost of folding. 
         */
        const std::string fold_formula;
        /// @brief The cost formula of the multicasting step for this layer.
        const std::string multicast_costs;
        /// @brief The src collapse formulation for the next layer.
        const std::string src_collapser;
        /// @brief The dst collapse formulation for the next layer.
        const std::string dst_collapser;
        /// @brief The context the layer is in.
        isl_ctx *const ctx;
}

int main(int argc, char* argv[])
{
    // Creates an isl context.
    isl_ctx *ctx = isl_ctx_alloc();

    // Creates the binding abstraction for the first layer.
    std::string srcs = R"SRC(
        {[id] -> [data] : id = 0 and data = id}
    )SRC";
    std::string data = R"DST(
        {[id, x, y] -> [data]: id = 0 and (-1 = x or x = 1) and 0 <= y <= 1 and data = y}
    )DST";
    binding test_case = binding(new binding_struct{srcs, data});

    // Calculates the cost formulas of the first layer.
    /// @note Read right to left like function composition.
    std::string crease_costs = "{ [id, x, y] -> x: x >= 0; [id, x, y] -> -x: x < 0 }";
    std::string fold_formula = "{ [id, x, y] -> [id, y] }";
    std::string multicast_costs = "{ [id, y] -> y+1 }";

    // Calculates the collapse formulas of the first layer.
    std::string dst_collapse_formula = "{ [id] -> [id, x, y] }";
    std::string src_collapse_formula = "{ [id] -> [id] }";
    collapse collapse_formulas = collapse(new collapse_struct{src_collapse_formula, dst_collapse_formula});

    BranchTwig test = BranchTwig(crease_costs, fold_formula, multicast_costs, collapse_formulas, ctx);
    std::cout << "Evaluating..." << std::endl;
    binding collapsed = test.evaluate(test_case);
    // Prints out the collapsed binding abstraction for the next layer.
    std::cout << "Collapsed: " << collapsed->srcs << std::endl;
    std::cout << "Missing: " << collapsed->dsts << std::endl;
    std::cout << "Done." << std::endl;

    ///@note Frees ctx to check for memory leaks through ISL.
    isl_ctx_free(ctx);
}