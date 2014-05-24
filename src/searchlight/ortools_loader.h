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
 * @file ortools_loader.h
 *
 * This file contains the definition of the or-tools loader class for the
 * model. or-tools does not export it, but we need it to create duplicate
 * models for the validator.
 *
 * Since the template definition is given here, the function will be defined
 * HERE, not in the constraint solver library. This means any changes in
 * the definition of the CPModelLoader class in the library might cause
 * binary compatibility issues and nasty bugs. Should be fine for this project
 * though.
 *
 * @author Alexander Kalinin
 * @author or-tools Team
 */

#ifndef SEARCHLIGHT_ORTOOLS_LOADER_H_
#define SEARCHLIGHT_ORTOOLS_LOADER_H_

namespace operations_research {
class CPModelLoader {
 public:
  explicit CPModelLoader(Solver* const solver) : solver_(solver) {}
  ~CPModelLoader() {}

  Solver* solver() const { return solver_; }

  // Builds integer expression from proto and stores it. It returns
  // true upon success.
  bool BuildFromProto(const CPIntegerExpressionProto& proto);
  // Builds constraint from proto and returns it.
  Constraint* BuildFromProto(const CPConstraintProto& proto);
  // Builds interval variable from proto and stores it. It returns
  // true upon success.
  bool BuildFromProto(const CPIntervalVariableProto& proto);
  // Builds sequence variable from proto and stores it. It returns
  // true upon success.
  bool BuildFromProto(const CPSequenceVariableProto& proto);

  // Returns stored integer expression.
  IntExpr* IntegerExpression(int index) const;
  // Returns stored interval variable.
  IntervalVar* IntervalVariable(int index) const;

  bool ScanOneArgument(int type_index, const CPArgumentProto& arg_proto,
                       int64* to_fill);

  bool ScanOneArgument(int type_index, const CPArgumentProto& arg_proto,
                       IntExpr** to_fill);

  bool ScanOneArgument(int type_index, const CPArgumentProto& arg_proto,
                       std::vector<int64>* to_fill);

  bool ScanOneArgument(int type_index, const CPArgumentProto& arg_proto,
                       IntTupleSet* to_fill);

  bool ScanOneArgument(int type_index, const CPArgumentProto& arg_proto,
                       std::vector<IntVar*>* to_fill);

  bool ScanOneArgument(int type_index, const CPArgumentProto& arg_proto,
                       IntervalVar** to_fill);

  bool ScanOneArgument(int type_index, const CPArgumentProto& arg_proto,
                       std::vector<IntervalVar*>* to_fill);

  bool ScanOneArgument(int type_index, const CPArgumentProto& arg_proto,
                       SequenceVar** to_fill);

  bool ScanOneArgument(int type_index, const CPArgumentProto& arg_proto,
                       std::vector<SequenceVar*>* to_fill);

  template <class P, class A>
  bool ScanArguments(const std::string& type, const P& proto, A* to_fill) {
    const int index = tags_.Index(type);
    for (int i = 0; i < proto.arguments_size(); ++i) {
      if (ScanOneArgument(index, proto.arguments(i), to_fill)) {
        return true;
      }
    }
    return false;
  }

  int TagIndex(const std::string& tag) const { return tags_.Index(tag); }

  void AddTag(const std::string& tag) { tags_.Add(tag); }

  // TODO(user): Use.
  void SetSequenceVariable(int index, SequenceVar* const var) {}

 private:
  Solver* const solver_;
  std::vector<IntExpr*> expressions_;
  std::vector<IntervalVar*> intervals_;
  std::vector<SequenceVar*> sequences_;
  VectorMap<std::string> tags_;
};
}
#endif /* SEARCHLIGHT_ORTOOLS_LOADER_H_ */
