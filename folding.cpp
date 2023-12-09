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

#include <barvinok/polylib.h>
#include <barvinok/barvinok.h>
#include <barvinok/isl.h>

/// @brief Strings representing the src and dst datum holds/requests in ISL.
struct geometry
{
    std::string srcs;
    std::string dsts;
};

class Layer
{
    public:
        /// @brief The cost formula of the folding step for this layer.
        const std::string crease_costs;
        /// @brief The folding action to a multicastable representation after
        /// calculating the cost of folding.
        const std::string fold_formula;
        /// @brief The cost formula of the multicasting step for this layer.
        const std::string multicast_costs;
        /// @brief The collapse formulation for the next layer.
        const std::string collapse_formula;
        /// @brief The context the layer is in.
        isl_ctx *const ctx;
    private:
        /// @brief This is the result of the cost for computations of this layer
        /// specifically, excluding the cost of the atomic units of this layer.
        double cost_result = 0;
        /// @brief Stores the cost of each operation in the layer.
        
        /// @brief This is the result of folding the cost of the atomic units
        /// of this layer into the folding formulation for the next layer.
        double latency = 0;
        std::string fold_result;
        /// @brief Whether this layer is dead.
        bool dead = false;
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
        Layer(
            const std::string& crease_costs, const std::string& fold_formula,
            const std::string& multicast_costs,
        isl_ctx* ctx):
        crease_costs(crease_costs), fold_formula(fold_formula),
        multicast_costs(multicast_costs), ctx(ctx) {}

        /// @brief Calculates the cost of the atomic units of this layer, then
        /// stores the cost.
        void evaluate(const std::string& s_srcs, const std::string& s_dsts)
        {
            // Folds the destinations onto their connected trunk.
            isl_map *folded = this->fold(s_dsts);
            std::cout << "Crease Cost: " << this->cost_result << std::endl;
            std::cout << "Folded: " << isl_map_to_str(folded) << std::endl;

            // Calculates the cost to every folded node per datum.
            this->multicast(folded);
            std::cout << "Total Cost: " << this->cost_result << std::endl;

            // Frees the maps.
        }
    private:
        /** @brief Folds the destinations onto their connected trunk. Updates
         * this->cost_result by adding the cost of casting to the folded dsts.
         * 
         * @param s_dsts The destinations to fold as an ISL string.
         * @return The folded destinations as an ISL string.
         * @post this->cost_result is incremented by the folding cost.
         */
        isl_map *fold(const std::string& dsts)
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
            // Adds the cost to the total cost.
            this->cost_result += isl_val_get_d(v_total_cost);
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

            return p_folded_condensed;
        }

        /** @brief Calculates the the cost to every folded node per datum and
         * adds it to this->cost_result.
         * 
         * @param folded_geometry The folded destinations from this->fold.
         * @post this->cost_result is incremented by the cost of multicasting.
         */
        void multicast(__isl_take isl_map *folded_geometry)
        {
            // Reads the multicast cost formulation into isl format.
            isl_pw_qpolynomial *p_cast_cost = isl_pw_qpolynomial_read_from_str(ctx, this->multicast_costs.c_str());
            // Applies the cost formulation to the folded dsts.
            isl_pw_qpolynomial *p_cost_applied = isl_map_apply_pw_qpolynomial(folded_geometry, p_cast_cost);
            // Sums all the costs.
            isl_pw_qpolynomial *p_total_cost = isl_pw_qpolynomial_sum(p_cost_applied);
            // Evaluates the cost.
            isl_val *v_total_cost = isl_pw_qpolynomial_eval(p_total_cost, isl_point_zero(isl_pw_qpolynomial_get_domain_space(p_total_cost)));
            // Adds the cost to the total cost.
            this->cost_result += isl_val_get_d(v_total_cost);
            // Frees the values.
            isl_val_free(v_total_cost);
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
    // isl_map *dst = isl_map_read_from_str(ctx, data.c_str());
    // std::string fold = R"FOLD(
    //     {[data] -> [id, y] : y = y'}
    // )FOLD";
    // std::string fold_data = "{ [id, y] -> [id, x, y] }";
    // std::string all_after = "{ [id, y] -> [id, y'] : y' > y }";

    // isl_map* id_to_all_x_y = isl_map_read_from_str(
    //     ctx,
    //     "{ [id] -> [id, x, y] }"
    // );
    // isl_map* id_to_all_dst_data = isl_map_apply_range(id_to_all_x_y, dst);
    // isl_map* id_to_missing_data = isl_map_subtract(id_to_all_dst_data, src);

    /// @todo Do abs not sq.
    std::string crease_costs = "{ [id, x, y] -> x^2 }";
    std::string fold_formula = "{ [id, x, y] -> [id, y] }";
    std::string multicast_costs = "{ [id, y] -> y }";
    Layer test = Layer(crease_costs, fold_formula, multicast_costs, ctx);
    std::cout << "Evaluating..." << std::endl;
    test.evaluate(srcs, data);

    ///@note Frees ctx to check for memory leaks through ISL.
    isl_ctx_free(ctx);
}