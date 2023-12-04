/** @brief This class represents a "layer" of analysis. It calculates the
 *  cost of the atomic units of this layer, then folds them into a problem
 * formulation for the next layer **/
#include "folding.h"
#include <memory>
#include <string>
#include <isl/aff.h>
#include <isl/map.h>

/// @brief Strings representing the src and dst datum holds/requests in ISL.
struct geometry
{
    std::string srcs;
    std::string dsts;
}

class Layer
{
    public:
        /// @brief The cost of the atomic units of this layer.
        const std::string cost_formula;
        /// @brief The folding formulation for the next layer.
        const std::string folding_formula;
        /// @brief The collapse formulation for the next layer.
        const std::string collapse_formula;
        /// @brief The context the layer is in.
        isl_ctx *const ctx;
    private:
        /// @brief This is the result of the cost for computations of this layer
        /// specifically, excluding the cost of the atomic units of this layer.
        const double cost_result = 0;
        /// @brief This is the result of folding the cost of the atomic units
        /// of this layer into the folding formulation for the next layer.
        double latency = 0;
        std::string fold_result;
        /// @brief Whether this layer is dead.
        bool dead = false;
    public:
        /// @brief Constructs a layer from a cost formulation and a folding
        /// formulation.
        Layer(std::string& cost_formula, std::string& folding_formula, isl_ctx* ctx):
        cost_formula(cost_formula), folding_formula(folding_formula), ctx(ctx) {}
        /// @brief Calculates the cost of the atomic units of this layer, then
        /// stores the cost.
        void evaluate(std::string s_srcs, std::string& s_dsts)
        {
            // Folds the destinations onto their connected trunk.
            isl_map *folded = this->fold(s_dsts);
            /** @todo Sum the cost of the atomic units of this layer and add to
              * the cost_result.*/

            // Calculates the cost to every folded node per datum.
            isl_map *expense = this->expense(folded);
            /** @todo Sum the cost of the expenses to each node and add to the
              * cost_result.*/

            // Formulates the result for the next layer.
            isl_map *collapsed = this->collapse(expense);
            // Reads the collapsed map into a string.
            fold_result = isl_map_to_str(collapsed);
            // Frees the maps.
            isl_map_free(folded);
            isl_map_free(expense);
            isl_map_free(collapsed);
        }
    private:
        /// @brief Folds the destinations onto their connected trunk.
        isl_map *fold(std::string& s_dsts)
        {
            // Reads the dsts into isl format.
            isl_map *dsts = isl_map_read_from_str(ctx, s_dsts.c_str());
            // Uses the folding formulation to calculate the cost of the units
            // in this layer to the bus trunk.
            isl_map *fold = isl_map_read_from_str(ctx, folding_formula.c_str());
            // Applies the folding formulation to the dsts.
            isl_map *dsts_folded = isl_map_apply_range(dsts, fold);
            
            return dsts_folded;
        }
        /// @brief Calculates the the cost to every folded node per datum.
        isl_map *expense(isl_map *fold)
        {
            // Reads the cost formulation into isl format.
            isl_map *cost = isl_map_read_from_str(ctx, cost_formula.c_str());
            // Applies the cost formulation to the folded dsts.
            isl_map *cost_applied = isl_map_apply_range(fold, cost);
            
            return cost_applied;
        }
};