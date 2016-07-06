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
#include <searchlight/relax.h>

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
void MimicAvg(Searchlight *sl, uint32_t id) {
    // SL Solver
    SearchlightSolver &sl_solver = sl->GetSLSolver(id);
    // CP solver
    Solver &solver = sl_solver.GetSearchSolver();

    // or-tools is deterministic by default
    solver.ReSeed(ACMRandom::HostnamePidTimeSeed());

    // parse the params
    const SearchlightConfig &config = sl->GetConfig();

    // task params
    const int32 start_id = config.get<int32>("mimic.l_id");
    const int32 end_id   = config.get<int32>("mimic.u_id");
    const int32 start_time = config.get<int32>("mimic.l_time");
    const int32 end_time   = config.get<int32>("mimic.u_time");
    const int32 avg_l   = config.get<int32>("mimic.avg_l");
    const int32 avg_u   = config.get<int32>("mimic.avg_u");
    const int32 avg_relax_l   = config.get<int32>("mimic.avg_relax_l", avg_l);
    const int32 avg_relax_u   = config.get<int32>("mimic.avg_relax_u", avg_u);
    const int32 len_l  = config.get<int32>("mimic.len_l");
    const int32 len_u  = config.get<int32>("mimic.len_u");
    const int32 step_len  = config.get<int32>("mimic.step_len", 1);
    const int32 step_time  = config.get<int32>("mimic.step_time", 1);

    const int32 time_limit = config.get("mimic.time_limit", 3600);
    const int luby_scale = config.get("searchlight.sl.luby_scale", 1);

    // problem params
    std::vector<IntVar *> coords(2);
    std::vector<IntVar *> lens(2);

    // contraction info
    SizeVector contr_constraints;
    BoolVector contr_spec;

    // coords
    coords[0] = MakeIntVarWithName(solver, start_id, end_id, 1, "id");
    coords[1] = MakeIntVarWithName(solver, start_time, end_time, step_time,
            "time");

    // lens
    lens[0] = solver.MakeIntConst(1);
    lens[1] = MakeIntVarWithName(solver, len_l, len_u, step_len, "len");

    // convenience -- all vars in a single vector
    std::vector<IntVar *> all_vars(coords);
    all_vars.insert(all_vars.end(), lens.begin(), lens.end());
    // vars expressed as IntExpr (required for UDFs)
    std::vector<IntExpr *> all_exprs(all_vars.begin(), all_vars.end());

    // Search vars
    std::vector<IntVar *> search_vars{coords[0], coords[1], lens[1]};

    // valid window
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords[1], lens[1]), end_time + 1));

    // average
    const std::string signal_name = config.get<std::string>("mimic.signal");
    AttributeID attr = sl->RegisterAttribute(signal_name);
    AdapterPtr adapter = sl_solver.CreateAdapter("mimic");
    UDFFunctionCreator avg_fab = sl->GetUDFFunctionCreator("avg");

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
        const int c_spec = config.get("mimic.avg_contr", 0);
        if (c_spec != 0) {
            contr_constraints.push_back(cid);
            contr_spec.push_back(c_spec == 1 ? true : false);
        }
        solver.AddConstraint(avg_bt);
    } else {
        solver.AddConstraint(solver.MakeBetweenCt(avg, avg_l, avg_u));
    }

    // neighborhood
    const int32 lneighb_size  = config.get("mimic.neighborhood.l_size", 0);
    const int32 rneighb_size  = config.get("mimic.neighborhood.r_size", 0);
    if (lneighb_size || rneighb_size) {
        /*
         * FIXME: This is a mess! And should be rewritten ASAP!
         */
        UDFFunctionCreator max_fab = sl->GetUDFFunctionCreator("max");
        std::vector<IntExpr *> parts(4);

        // max of the region
        IntExpr * const region_max = solver.RevAlloc(max_fab(&solver, adapter,
                all_exprs, udf_params));

        // left
        if (lneighb_size) {
            parts[0] = coords[0]; // the same id/timeline
            parts[1] = solver.MakeSum(coords[1], -lneighb_size);
            parts[2] = solver.MakeIntConst(1);
            parts[3] = solver.MakeIntConst(lneighb_size);

            // valid left neighborhood
            solver.AddConstraint(solver.MakeGreaterOrEqual(parts[1],
                    start_time));

            // difference of maximums
            const int32 left_max_diff  = config.get(
                    "mimic.neighborhood.left_max_diff", 0);
            const int32 left_relax_diff  = config.get(
                    "mimic.neighborhood.left_relax_diff", left_max_diff);
            // expressions
            IntVar * const left_max = solver.RevAlloc(max_fab(&solver,
                    adapter, parts, udf_params))->Var();
            IntExpr * const rl_max_diff = solver.MakeDifference(
                    region_max, left_max);
            sl->AddTrackExpr(rl_max_diff, "left_diff");
            if (config.get("relax.on", false)) {
                RelaxableConstraint *l_diff_gt = solver.RevAlloc(
                        new searchlight::GreaterEqExprCst(
                            &solver, rl_max_diff, left_max_diff));
                const size_t cid = sl->RegisterConstraint("left_diff_const", id,
                                                          l_diff_gt,
                                                          left_relax_diff,
                                                          left_relax_diff);
                const int c_spec = config.get("mimic.left_contr", 0);
                if (c_spec != 0) {
                    contr_constraints.push_back(cid);
                    contr_spec.push_back(c_spec == 1 ? true : false);
                }
                solver.AddConstraint(l_diff_gt);
            } else {
                solver.AddConstraint(solver.MakeGreater(rl_max_diff,
                        left_max_diff));
            }
        }

        // right
        if (rneighb_size) {
            // right
            parts[0] = coords[0]; // the same id/timeline
            parts[1] = solver.MakeSum(coords[1], lens[1]);
            parts[2] = solver.MakeIntConst(1);
            parts[3] = solver.MakeIntConst(rneighb_size);

            // valid right neighborhood
            solver.AddConstraint(solver.MakeLessOrEqual(
                    solver.MakeSum(parts[1], rneighb_size),
                    end_time + 1));

            // difference of maximums
            const int32 right_max_diff = config.get(
                    "mimic.neighborhood.right_max_diff", 0);
            const int32 right_relax_diff  = config.get(
                    "mimic.neighborhood.right_relax_diff", right_max_diff);
            // expressions
            IntVar * const right_max = solver.RevAlloc(max_fab(&solver,
                    adapter, parts, udf_params))->Var();
            IntExpr * const rr_max_diff = solver.MakeDifference(
                    region_max, right_max);
            sl->AddTrackExpr(rr_max_diff, "right_diff");
            if (config.get("relax.on", false)) {
                RelaxableConstraint *r_diff_gt = solver.RevAlloc(
                        new searchlight::GreaterEqExprCst(
                            &solver, rr_max_diff, right_max_diff));
                const size_t cid = sl->RegisterConstraint("right_diff_const",
                                                          id, r_diff_gt,
                                                          right_relax_diff,
                                                          right_relax_diff);
                const int c_spec = config.get("mimic.right_contr", 0);
                if (c_spec != 0) {
                    contr_constraints.push_back(cid);
                    contr_spec.push_back(c_spec == 1 ? true : false);
                }
                solver.AddConstraint(r_diff_gt);
            } else {
                solver.AddConstraint(solver.MakeGreater(rr_max_diff,
                        right_max_diff));
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
    const std::string search_heuristic = config.get<std::string>("mimic.db");
    DecisionBuilder *db;
    std::vector<SearchMonitor *> mons;
    if (search_heuristic == "impact") {
        db = solver.MakeDefaultPhase(search_vars);
    } else if (search_heuristic == "random") {
        db = solver.MakePhase(search_vars, Solver::CHOOSE_RANDOM,
            Solver::ASSIGN_RANDOM_VALUE);
        if (luby_scale != 0) {
            mons.push_back(solver.MakeLubyRestart(luby_scale));
        }
        mons.push_back(solver.MakeTimeLimit(time_limit * 1000));
    } else if (search_heuristic == "split") {
        db = solver.MakePhase(search_vars, Solver::CHOOSE_MAX_SIZE,
            Solver::SPLIT_LOWER_HALF);
        const bool balance = config.get("balance.solver_balance", 1);
        const double low_thr = config.get("balance.general_low", 0.1);
        const double high_thr = config.get("balance.general_high", 0.5);
        if (balance) {
            mons.push_back(sl_solver.CreateBalancingMonitor(search_vars,
                    low_thr, high_thr));
        }
    } else if (search_heuristic == "sl") {
        if (time_limit != 0) {
            mons.push_back(MakeCumulativeTimeLimit(solver, time_limit * 1000));
        }
        db = sl_solver.CreateDefaultHeuristic(coords, {lens[1]});
    } else {
        db = solver.MakePhase(search_vars, Solver::CHOOSE_FIRST_UNBOUND,
            Solver::ASSIGN_MIN_VALUE);
    }

    // set task
    sl_solver.SetTask(coords, {lens[1]}, db, mons);
}
} /* namespace searchlight */
