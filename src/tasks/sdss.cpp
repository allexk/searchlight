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
 * @file sdss.cpp
 *
 * Implementation of the SDSS task.
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
IntVar *MakeIntVarWithName(Solver &solver, int64 start, int64 end,
        int64 step, const std::string &name) {
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

std::vector<double> ReadJSONArray(
        const searchlight::SearchlightConfig &sl_config,
        const std::string &path) {
    std::vector<double> res;
    try {
        for (const auto &node: sl_config.get_child(path)) {
            const int val = node.second.get_value<double>();
            res.push_back(val);
        }
    } catch (const boost::property_tree::ptree_bad_path &) {
        std::cerr << "Average values for ugriz are not specified!";
        throw;
    }
    return res;
}

}

extern "C"
void SdssUgrizAvg(Searchlight *sl) {
    // get the solver
    Solver &solver = sl->GetSolver();

    // or-tools is deterministic by default
    solver.ReSeed(ACMRandom::HostnamePidTimeSeed());

    // parse the params
    const SearchlightConfig &config = sl->GetConfig();

    // task params
    const int64 start_x = config.get<int64>("sdss.lx");
    const int64 end_x   = config.get<int64>("sdss.ux");
    const int64 start_y = config.get<int64>("sdss.ly");
    const int64 end_y   = config.get<int64>("sdss.uy");
    const int64 size_l  = config.get<int64>("sdss.size_l");
    const int64 size_u  = config.get<int64>("sdss.size_u");
    const int64 len_lx  = config.get<int64>("sdss.len_lx");
    const int64 len_ux  = config.get<int64>("sdss.len_ux");
    const int64 len_ly  = config.get<int64>("sdss.len_ly");
    const int64 len_uy  = config.get<int64>("sdss.len_uy");
    const int64 step_x  = config.get<int64>("sdss.step_x", 1);
    const int64 step_y  = config.get<int64>("sdss.step_y", 1);
    const int64 step_lx  = config.get<int64>("sdss.step_lx", 1);
    const int64 step_ly  = config.get<int64>("sdss.step_ly", 1);

    auto avg_low = ReadJSONArray(config, "sdss.avg_l");
    auto avg_high = ReadJSONArray(config, "sdss.avg_h");

    auto min_low = ReadJSONArray(config, "sdss.min_l");
    auto min_high = ReadJSONArray(config, "sdss.min_h");

    const int32 time_limit = config.get("sdss.time_limit", 3600);
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

    // valid window
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords[0], lens[0]), end_x + 1));
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords[1], lens[1]), end_y + 1));
    solver.AddConstraint(solver.MakeBetweenCt(
            solver.MakeProd(lens[0], lens[1]), size_l, size_u));

    // average: we have 5 attributes: u,g,r,i,z
    const char *ugriz_attrs[]{"u", "g", "r", "i", "z"};
    AdapterPtr adapter = sl->CreateAdapter(std::string("sw"));
    for (int i = 0; i < 5; i++) {
        const double l = avg_low[i];
        const double h = avg_high[i];
        const double ml = min_low[i];
        const double mh = min_high[i];
        if (l <= h || ml <= mh) {
            const std::string attr_str = std::string(ugriz_attrs[i]);
            const AttributeID attr = sl->RegisterAttribute(attr_str);

            // avg
            UDFFunctionCreator avg_fab = sl->GetUDFFunctionCreator("avg");
            std::vector<int64> udf_params(1, int64(attr));
            IntExpr * const avg = solver.RevAlloc(avg_fab(
                    &solver, adapter, all_vars, udf_params));

            solver.AddConstraint(solver.MakeBetweenCt(avg, l, h));

            // min
            UDFFunctionCreator min_fab = sl->GetUDFFunctionCreator("min");
            IntExpr * const min = solver.RevAlloc(min_fab(
                    &solver, adapter, all_vars, udf_params));

            solver.AddConstraint(solver.MakeGreaterOrEqual(min, int64(ml)));

            // max
            UDFFunctionCreator max_fab = sl->GetUDFFunctionCreator("max");
            IntExpr * const max = solver.RevAlloc(max_fab(
                    &solver, adapter, all_vars, udf_params));

            solver.AddConstraint(solver.MakeLessOrEqual(max, int64(mh)));
        }


    }

    // create the search phase
    const std::string search_heuristic = config.get<std::string>("sdss.db");
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
    } else if (search_heuristic == "sl") {
        if (time_limit != 0) {
            mons.push_back(MakeCumulativeTimeLimit(solver, time_limit * 1000));
        }
        db = sl->CreateDefaultHeuristic(coords, lens);
    } else {
        db = solver.MakePhase(all_vars, Solver::CHOOSE_FIRST_UNBOUND,
            Solver::ASSIGN_MIN_VALUE);
    }

    // solve!
    sl->Prepare(coords, lens, db, mons);
}
} /* namespace searchlight */
