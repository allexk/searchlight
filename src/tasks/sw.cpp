/* Copyright 2014, Brown University, Providence, RI.
 *
 *                         All Rights Reserved
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose other than its incorporation into a
 * commercial product is hereby granted without fee, provided that the
 * above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation, and that the name of Brown University not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific, written prior permission.
 *
 * BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
 * INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
 * PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
 * ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

/**
 * @file sw.cpp
 *
 * Implementation of Semantic Windows (SW) task.
 *
 * @author Alexander Kalinin
 */

#include <searchlight/searchlight.h>

#include <boost/lexical_cast.hpp>

namespace searchlight {

namespace {

/**
 * Creates an integer variable and assigns it a name.
 *
 * The deal here is that the solver may optimize the resulting variable by,
 * for example, altering its domain. It will still return a "correct" variable,
 * but the variable might be a cast variable for some internal expression.
 * Since names are really important in Searchlight for identification, we
 * enforce the name after the variable is created.
 *
 * @param solver solver to create the variable for
 * @param start left boundary for the domain (inclusive)
 * @param end right boundary for the domain (inclusive)
 * @param step step for values in the domain (1 means all values)
 * @param name the name for the variable
 * @return integer variable with the specified name
 */
IntVar *MakeIntVarWithName(Solver &solver, int32 start, int32 end,
        int32 step, const std::string &name) {
    IntVar *var;
    if (step == 1) {
        var = solver.MakeIntVar(start, end);
    } else {
        std::vector<int64> vals;
        for (int64 v = start; v <= end; v += step) {
            vals.push_back(v);
        }
        var = solver.MakeIntVar(vals);
    }
    var->set_name(name);
    return var;
}

}

extern "C"
void SemWindowsAvg(Searchlight *sl, uint32_t id) {
    // SL Solver
    SearchlightSolver &sl_solver = sl->GetSLSolver(id);
    // CP solver
    Solver &solver = sl_solver.GetSearchSolver();

    // or-tools is deterministic by default
    solver.ReSeed(ACMRandom::HostnamePidTimeSeed());

    // parse the params
    const SearchlightConfig &config = sl->GetConfig();

    // task params
    const int32 start_x = config.get<int32>("sw.lx");
    const int32 end_x   = config.get<int32>("sw.ux");
    const int32 start_y = config.get<int32>("sw.ly");
    const int32 end_y   = config.get<int32>("sw.uy");
    const int32 avg_l   = config.get<int32>("sw.avg_l");
    const int32 avg_u   = config.get<int32>("sw.avg_u");
    const int32 avg_relax_l   = config.get<int32>("sw.avg_relax_l", avg_l);
    const int32 avg_relax_u   = config.get<int32>("sw.avg_relax_u", avg_u);
    const int32 size_l  = config.get<int32>("sw.size_l");
    const int32 size_u  = config.get<int32>("sw.size_u");
    const int32 len_lx  = config.get<int32>("sw.len_lx");
    const int32 len_ux  = config.get<int32>("sw.len_ux");
    const int32 len_ly  = config.get<int32>("sw.len_ly");
    const int32 len_uy  = config.get<int32>("sw.len_uy");
    const int32 step_x  = config.get<int32>("sw.step_x", 1);
    const int32 step_y  = config.get<int32>("sw.step_y", 1);
    const int32 step_lx  = config.get<int32>("sw.step_lx", 1);
    const int32 step_ly  = config.get<int32>("sw.step_ly", 1);

    const int32 time_limit = config.get("sw.time_limit", 3600);
    const int luby_scale = config.get("searchlight.sl.luby_scale", 1);

    // problem params
    std::vector<IntVar *> coords(2);
    std::vector<IntVar *> lens(2);

    // coords
    coords[0] = MakeIntVarWithName(solver, start_x, end_x, step_x, "x");
    coords[1] = MakeIntVarWithName(solver, start_y, end_y, step_y, "y");

    // lens
    lens[0] = MakeIntVarWithName(solver, len_lx, len_ux, step_lx, "lx");
    lens[1] = MakeIntVarWithName(solver, len_ly, len_uy, step_ly, "ly");

    // convenience -- all vars in a single vector
    std::vector<IntVar *> all_vars(coords);
    all_vars.insert(all_vars.end(), lens.begin(), lens.end());
    // vars expressed as IntExpr (required for UDFs)
    std::vector<IntExpr *> all_exprs(all_vars.begin(), all_vars.end());

    // valid window
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords[0], lens[0]), end_x + 1));
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords[1], lens[1]), end_y + 1));
    solver.AddConstraint(solver.MakeBetweenCt(
            solver.MakeProd(lens[0], lens[1]), size_l, size_u));

    // average
    AttributeID attr = sl->RegisterAttribute("val");
    AdapterPtr adapter = sl_solver.CreateAdapter("sw");
    UDFFunctionCreator avg_fab = sl->GetUDFFunctionCreator("avg");

    // contraction info
    SizeVector contr_constraints;
    BoolVector contr_spec;

    std::vector<int64> udf_params(1, int64(attr));
    IntExpr * const avg = solver.RevAlloc(avg_fab(&solver, adapter, all_exprs,
            udf_params));
    sl->AddTrackExpr(avg, "avg");
    if (config.get("relax.on", false)) {
        RelaxableConstraint *avg_bt = solver.RevAlloc(
                new searchlight::BetweenCt(
                    &solver, avg, avg_l, avg_u));
        const size_t cid = sl->RegisterConstraint("avg_const", id, avg_bt,
                                                  avg_relax_l, avg_relax_u);
        const int c_spec = config.get("sw.avg_contr", 0);
        if (c_spec != 0) {
            contr_constraints.push_back(cid);
            contr_spec.push_back(c_spec == 1 ? true : false);
        }
        solver.AddConstraint(avg_bt);
    } else {
        solver.AddConstraint(solver.MakeBetweenCt(avg, avg_l, avg_u));
    }

    // neighborhood
    const int32 neighb_size  = config.get("sw.neighborhood.size", 0);
    if (neighb_size) {
        /*
         * FIXME: This is a mess! And should be rewritten ASAP!
         */
        UDFFunctionCreator max_fab = sl->GetUDFFunctionCreator("max");
        UDFFunctionCreator min_fab = sl->GetUDFFunctionCreator("min");
        std::vector<IntExpr *> parts(4);
        std::vector<IntExpr *> part_maxs(4);

        // left
        parts[0] = solver.MakeSum(coords[0], -neighb_size);
        parts[1] = coords[1];
        parts[2] = solver.MakeIntConst(neighb_size);
        parts[3] = lens[1];
        part_maxs[0] = solver.RevAlloc(max_fab(&solver,
                adapter, parts, udf_params));

        // top
        parts[0] = solver.MakeSum(coords[0], -neighb_size);
        parts[1] = solver.MakeSum(coords[1], -neighb_size);
        parts[2] = solver.MakeSum(lens[0], 2 * neighb_size);
        parts[3] = solver.MakeIntConst(neighb_size);
        part_maxs[1] = solver.RevAlloc(max_fab(&solver,
                adapter, parts, udf_params));

        // right
        parts[0] = solver.MakeSum(coords[0], lens[0]);
        parts[1] = coords[1];
        parts[2] = solver.MakeIntConst(neighb_size);
        parts[3] = lens[1];
        part_maxs[2] = solver.RevAlloc(max_fab(&solver,
                adapter, parts, udf_params));

        // bottom
        parts[0] = solver.MakeSum(coords[0], -neighb_size);
        parts[1] = solver.MakeSum(coords[1], lens[1]);
        parts[2] = solver.MakeSum(lens[0], 2 * neighb_size);
        parts[3] = solver.MakeIntConst(neighb_size);
        part_maxs[3] = solver.RevAlloc(max_fab(&solver,
                adapter, parts, udf_params));

        /*
         * Max of the neighborhood. We do the explicit nested maxs here to
         * avoid cast constraints. Cast constraints mess up the relaxation
         * process.
         */
        IntExpr * const nmax = solver.MakeMax(
                solver.MakeMax(part_maxs[0], part_maxs[1]),
                solver.MakeMax(part_maxs[2], part_maxs[3]));

        const int32 max_diff  = config.get("sw.neighborhood.max_diff", 0);
        const int32 max_relax_diff  = config.get(
                "sw.neighborhood.max_relax_diff", max_diff);
        const int32 max_diff_high  = config.get(
                "sw.neighborhood.max_diff_high", max_diff + 1000);
        if (max_diff) {
            // max of the region
            IntExpr * const rmax = solver.RevAlloc(max_fab(&solver, adapter,
                    all_exprs, udf_params));
            IntExpr * const rn_max_diff = solver.MakeDifference(rmax, nmax);
            sl->AddTrackExpr(rn_max_diff, "max_diff");
            if (config.get("relax.on", false)) {
                RelaxableConstraint *max_diff_gt = solver.RevAlloc(
                        new searchlight::GreaterEqExprCst(
                            &solver, rn_max_diff, max_diff));
                const size_t cid = sl->RegisterConstraint("max_diff_const", id,
                                                          max_diff_gt,
                                                          max_relax_diff,
                                                          max_diff_high);
                const int c_spec = config.get("sw.neighborhood.max_diff_contr",
                                              0);
                if (c_spec != 0) {
                    contr_constraints.push_back(cid);
                    contr_spec.push_back(c_spec == 1 ? true : false);
                }
                solver.AddConstraint(max_diff_gt);
            } else {
                solver.AddConstraint(solver.MakeGreater(rn_max_diff, max_diff));
            }
        }

        const int32 min_diff  = config.get("sw.neighborhood.min_diff", 0);
        const int32 min_relax_diff  = config.get(
                "sw.neighborhood.min_relax_diff", min_diff);
        const int32 min_diff_high  = config.get(
                "sw.neighborhood.min_diff_high", min_diff + 1000);
        if (min_diff) {
            // min of the region
            IntExpr * const rmin = solver.RevAlloc(min_fab(&solver, adapter,
                    all_exprs, udf_params));
            IntExpr * const rn_min_diff = solver.MakeDifference(rmin, nmax);
            sl->AddTrackExpr(rn_min_diff, "min_diff");
            if (config.get("relax.on", false)) {
                RelaxableConstraint *min_diff_gt = solver.RevAlloc(
                        new searchlight::GreaterEqExprCst(
                            &solver, rn_min_diff, min_diff));
                const size_t cid = sl->RegisterConstraint("min_diff_const", id,
                                                          min_diff_gt,
                                                          min_relax_diff,
                                                          min_diff_high);
                const int c_spec = config.get("sw.neighborhood.min_diff_contr",
                                              0);
                if (c_spec != 0) {
                    contr_constraints.push_back(cid);
                    contr_spec.push_back(c_spec == 1 ? true : false);
                }
                solver.AddConstraint(min_diff_gt);
            } else {
                solver.AddConstraint(solver.MakeGreater(rn_min_diff, min_diff));
            }
        }
    }

    // Enable contraction
    const int contr_type = config.get("relax.contr_type", 0);
    if (!contr_constraints.empty() && contr_type) {
        Relaxator::ContractionType type;
        type = contr_type == 2 ? Relaxator::ContractionType::SKYLINE :
                Relaxator::ContractionType::RANK;
        sl_solver.EnableContraction(contr_constraints, contr_spec, type);
    }

    // create the search phase
    const std::string search_heuristic = config.get<std::string>("sw.db");
    DecisionBuilder *db;
    std::vector<SearchMonitor *> mons;
    if (search_heuristic == "impact") {
        db = solver.MakeDefaultPhase(all_vars);
    } else if (search_heuristic == "random") {
        db = solver.MakePhase(all_vars, Solver::CHOOSE_RANDOM,
            Solver::ASSIGN_RANDOM_VALUE);
        if (luby_scale != 0) {
            mons.push_back(solver.MakeLubyRestart(luby_scale));
        }
        mons.push_back(solver.MakeTimeLimit(time_limit * 1000));
    } else if (search_heuristic == "split") {
        db = solver.MakePhase(all_vars, Solver::CHOOSE_MAX_SIZE,
            Solver::SPLIT_LOWER_HALF);
        const bool balance = config.get("balance.solver_balance", 1);
        const double low_thr = config.get("balance.general_low", 0.1);
        const double high_thr = config.get("balance.general_high", 0.5);
        if (balance) {
            mons.push_back(sl_solver.CreateBalancingMonitor(all_vars, low_thr,
                    high_thr));
        }
    } else if (search_heuristic == "sl") {
        if (time_limit != 0) {
            mons.push_back(MakeCumulativeTimeLimit(solver, time_limit * 1000));
        }
        db = sl_solver.CreateDefaultHeuristic(coords, lens);
    } else {
        db = solver.MakePhase(all_vars, Solver::CHOOSE_FIRST_UNBOUND,
            Solver::ASSIGN_MIN_VALUE);
    }

    // set task
    sl_solver.SetTask(coords, lens, db, mons);
}
} /* namespace searchlight */
