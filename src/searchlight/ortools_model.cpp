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
 * @file ortools_model.cpp
 *
 * This file contains copy-pasted and slightly modified original or-tools
 * classes to manipulate CP models, for example to clone them.
 *
 * Most (if not all) of this stuff is already present in the or-tools
 * implementation, but it is not exported for general use. Plus, some of the
 * facilities miss important properties, like packing all variable names
 * during model export. For now, only real variable names are packed, but not
 * optimized variables, like x + c or x * c, which happen quite frequently.
 *
 * Comments for every copy-pasted class are given before the class definition.
 *
 * @author Alexander Kalinin
 * @author or-tools Team
 */

#include "ortools_model.h"
#include "relax.h"

#include <base/stl_util.h>
#include <base/callback.h>

/*
 * Copy-pasted visitors for exporting models. We slightly modify SecondPass
 * visitor to export names for all, not only main, decision variables.
 */
namespace operations_research {
namespace {

// ---------- Model Protobuf Writers -----------

// ----- First Pass visitor -----

// This visitor collects all constraints and expressions.  It sorts the
// expressions, such that we can build them in sequence using
// previously created expressions.
class FirstPassVisitor : public ModelVisitor {
 public:
  FirstPassVisitor() {}  // Needed for Visual Studio.
  virtual ~FirstPassVisitor() {}

  virtual std::string DebugString() const { return "FirstPassVisitor"; }

  // Begin/End visit element.
  virtual void BeginVisitModel(const std::string& solver_name) {
    // Reset statistics.
    expression_map_.clear();
    delegate_map_.clear();
    expression_list_.clear();
    constraint_list_.clear();
    interval_list_.clear();
    sequence_list_.clear();
  }

  virtual void EndVisitConstraint(const std::string& type_name,
                                  const Constraint* const constraint) {
    Register(constraint);
  }

  virtual void EndVisitIntegerExpression(const std::string& type_name,
                                         const IntExpr* const expression) {
    Register(expression);
  }

  virtual void VisitIntegerVariable(const IntVar* const variable,
                                    IntExpr* const delegate) {
    if (delegate != nullptr) {
      delegate->Accept(this);
      delegate_map_[variable] = delegate;
    }
    Register(variable);
  }

  virtual void VisitIntegerVariable(const IntVar* const variable,
                                    const std::string& operation, int64 value,
                                    IntVar* const delegate) {
    delegate->Accept(this);
    delegate_map_[variable] = delegate;
    Register(variable);
  }

  virtual void VisitIntervalVariable(const IntervalVar* const variable,
                                     const std::string& operation, int64 value,
                                     IntervalVar* const delegate) {
    if (delegate != nullptr) {
      delegate->Accept(this);
    }
    Register(variable);
  }

  virtual void VisitSequenceVariable(const SequenceVar* const sequence) {
    for (int i = 0; i < sequence->size(); ++i) {
      sequence->Interval(i)->Accept(this);
    }
    Register(sequence);
  }

  // Visit integer expression argument.
  virtual void VisitIntegerExpressionArgument(const std::string& arg_name,
                                              IntExpr* const argument) {
    VisitSubArgument(argument);
  }

  virtual void VisitIntegerVariableArrayArgument(
      const std::string& arg_name, const std::vector<IntVar*>& arguments) {
    for (int i = 0; i < arguments.size(); ++i) {
      VisitSubArgument(arguments[i]);
    }
  }

  // Visit interval argument.
  virtual void VisitIntervalArgument(const std::string& arg_name,
                                     IntervalVar* const argument) {
    VisitSubArgument(argument);
  }

  virtual void VisitIntervalArrayArgument(
      const std::string& arg_name, const std::vector<IntervalVar*>& arguments) {
    for (int i = 0; i < arguments.size(); ++i) {
      VisitSubArgument(arguments[i]);
    }
  }

  // Visit sequence argument.
  virtual void VisitSequenceArgument(const std::string& arg_name,
                                     SequenceVar* const argument) {
    VisitSubArgument(argument);
  }

  virtual void VisitSequenceArrayArgument(
      const std::string& arg_name, const std::vector<SequenceVar*>& arguments) {
    for (int i = 0; i < arguments.size(); ++i) {
      VisitSubArgument(arguments[i]);
    }
  }

  // Export
  const hash_map<const IntExpr*, int>& expression_map() const {
    return expression_map_;
  }
  const hash_map<const IntervalVar*, int>& interval_map() const {
    return interval_map_;
  }
  const hash_map<const SequenceVar*, int>& sequence_map() const {
    return sequence_map_;
  }
  const hash_map<const IntVar*, IntExpr*>& delegate_map() const {
    return delegate_map_;
  }
  const std::vector<const IntExpr*>& expression_list() const {
    return expression_list_;
  }
  const std::vector<const Constraint*>& constraint_list() const {
    return constraint_list_;
  }
  const std::vector<const IntervalVar*>& interval_list() const {
    return interval_list_;
  }
  const std::vector<const SequenceVar*>& sequence_list() const {
    return sequence_list_;
  }

 private:
  void Register(const IntExpr* const expression) {
    if (!ContainsKey(expression_map_, expression)) {
      const int index = expression_map_.size();
      CHECK_EQ(index, expression_list_.size());
      expression_map_[expression] = index;
      expression_list_.push_back(expression);
    }
  }

  void Register(const Constraint* const constraint) {
    constraint_list_.push_back(constraint);
  }

  void Register(const IntervalVar* const interval) {
    if (!ContainsKey(interval_map_, interval)) {
      const int index = interval_map_.size();
      CHECK_EQ(index, interval_list_.size());
      interval_map_[interval] = index;
      interval_list_.push_back(interval);
    }
  }

  void Register(const SequenceVar* const sequence) {
    if (!ContainsKey(sequence_map_, sequence)) {
      const int index = sequence_map_.size();
      CHECK_EQ(index, sequence_list_.size());
      sequence_map_[sequence] = index;
      sequence_list_.push_back(sequence);
    }
  }

  void VisitSubArgument(IntExpr* const expression) {
    if (!ContainsKey(expression_map_, expression)) {
      expression->Accept(this);
    }
  }

  void VisitSubArgument(IntervalVar* const interval) {
    if (!ContainsKey(interval_map_, interval)) {
      interval->Accept(this);
    }
  }

  void VisitSubArgument(SequenceVar* const sequence) {
    if (!ContainsKey(sequence_map_, sequence)) {
      sequence->Accept(this);
    }
  }

  const std::string filename_;
  hash_map<const IntExpr*, int> expression_map_;
  hash_map<const IntervalVar*, int> interval_map_;
  hash_map<const SequenceVar*, int> sequence_map_;
  hash_map<const IntVar*, IntExpr*> delegate_map_;
  std::vector<const IntExpr*> expression_list_;
  std::vector<const Constraint*> constraint_list_;
  std::vector<const IntervalVar*> interval_list_;
  std::vector<const SequenceVar*> sequence_list_;
};

// ----- Argument Holder -----

class ArgumentHolder {
 public:
  template <class P>
  void ExportToProto(VectorMap<std::string>* const tags, P* const proto) const {
    for (const auto& it : integer_argument_) {
      CPArgumentProto* const arg_proto = proto->add_arguments();
      arg_proto->set_argument_index(tags->Add(it.first));
      arg_proto->set_integer_value(it.second);
    }

    for (const auto& it : integer_array_argument_) {
      CPArgumentProto* const arg_proto = proto->add_arguments();
      arg_proto->set_argument_index(tags->Add(it.first));
      for (int64 value : it.second) {
        arg_proto->add_integer_array(value);
      }
    }

    for (const auto& it : integer_matrix_argument_) {
      CPArgumentProto* const arg_proto = proto->add_arguments();
      arg_proto->set_argument_index(tags->Add(it.first));
      CPIntegerMatrixProto* const matrix_proto =
          arg_proto->mutable_integer_matrix();
      const int columns = it.second.first;
      CHECK_GT(columns, 0);
      const int rows = it.second.second.size() / columns;
      matrix_proto->set_rows(rows);
      matrix_proto->set_columns(columns);
      for (int64 value : it.second.second) {
        matrix_proto->add_values(value);
      }
    }

    for (const auto& it : integer_expression_argument_) {
      CPArgumentProto* const arg_proto = proto->add_arguments();
      arg_proto->set_argument_index(tags->Add(it.first));
      arg_proto->set_integer_expression_index(it.second);
    }

    for (const auto& it : integer_variable_array_argument_) {
      CPArgumentProto* const arg_proto = proto->add_arguments();
      arg_proto->set_argument_index(tags->Add(it.first));
      for (int expr : it.second) {
        arg_proto->add_integer_expression_array(expr);
      }
    }

    for (const auto& it : interval_argument_) {
      CPArgumentProto* const arg_proto = proto->add_arguments();
      arg_proto->set_argument_index(tags->Add(it.first));
      arg_proto->set_interval_index(it.second);
    }

    for (const auto& it : interval_array_argument_) {
      CPArgumentProto* const arg_proto = proto->add_arguments();
      arg_proto->set_argument_index(tags->Add(it.first));
      for (int arg : it.second) {
        arg_proto->add_interval_array(arg);
      }
    }

    for (const auto& it : sequence_argument_) {
      CPArgumentProto* const arg_proto = proto->add_arguments();
      arg_proto->set_argument_index(tags->Add(it.first));
      arg_proto->set_sequence_index(it.second);
    }

    for (const auto& it : sequence_array_argument_) {
      CPArgumentProto* const arg_proto = proto->add_arguments();
      arg_proto->set_argument_index(tags->Add(it.first));
      for (int arg : it.second) {
        arg_proto->add_sequence_array(arg);
      }
    }
  }

  const std::string& type_name() const { return type_name_; }

  void set_type_name(const std::string& type_name) { type_name_ = type_name; }

  void set_integer_argument(const std::string& arg_name, int64 value) {
    integer_argument_[arg_name] = value;
  }

  void set_integer_array_argument(const std::string& arg_name,
                                  const std::vector<int64>& values) {
    integer_array_argument_[arg_name] = values;
  }

  void set_integer_matrix_argument(const std::string& arg_name,
                                   const IntTupleSet& values) {
    const int rows = values.NumTuples();
    const int columns = values.Arity();
    std::pair<int, std::vector<int64>> matrix = std::make_pair(columns, std::vector<int64>());
    integer_matrix_argument_[arg_name] = matrix;
    std::vector<int64>* const vals = &integer_matrix_argument_[arg_name].second;
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < columns; ++j) {
        vals->push_back(values.Value(i, j));
      }
    }
  }

  void set_integer_expression_argument(const std::string& arg_name, int index) {
    integer_expression_argument_[arg_name] = index;
  }

  void set_integer_variable_array_argument(const std::string& arg_name,
                                           const int* const indices, int size) {
    for (int i = 0; i < size; ++i) {
      integer_variable_array_argument_[arg_name].push_back(indices[i]);
    }
  }

  void set_interval_argument(const std::string& arg_name, int index) {
    interval_argument_[arg_name] = index;
  }

  void set_interval_array_argument(const std::string& arg_name,
                                   const int* const indices, int size) {
    for (int i = 0; i < size; ++i) {
      interval_array_argument_[arg_name].push_back(indices[i]);
    }
  }

  void set_sequence_argument(const std::string& arg_name, int index) {
    sequence_argument_[arg_name] = index;
  }

  void set_sequence_array_argument(const std::string& arg_name,
                                   const int* const indices, int size) {
    for (int i = 0; i < size; ++i) {
      sequence_array_argument_[arg_name].push_back(indices[i]);
    }
  }

  int64 FindIntegerArgumentWithDefault(const std::string& arg_name, int64 def) {
    return FindWithDefault(integer_argument_, arg_name, def);
  }

  int64 FindIntegerArgumentOrDie(const std::string& arg_name) {
    return FindOrDie(integer_argument_, arg_name);
  }

  int64 FindIntegerExpressionArgumentOrDie(const std::string& arg_name) {
    return FindOrDie(integer_expression_argument_, arg_name);
  }

 private:
  std::string type_name_;
  hash_map<std::string, int> integer_expression_argument_;
  hash_map<std::string, int64> integer_argument_;
  hash_map<std::string, int> interval_argument_;
  hash_map<std::string, int> sequence_argument_;
  hash_map<std::string, std::vector<int64>> integer_array_argument_;
  hash_map<std::string, std::pair<int, std::vector<int64>>> integer_matrix_argument_;
  hash_map<std::string, std::vector<int>> integer_variable_array_argument_;
  hash_map<std::string, std::vector<int>> interval_array_argument_;
  hash_map<std::string, std::vector<int>> sequence_array_argument_;
};

// ----- Second Pass Visitor -----

static const int kModelVersion = 1;

// The second pass visitor will visited sorted expressions, interval
// vars and expressions and export them to a CPModelProto protocol
// buffer.
class SecondPassVisitor : public ModelVisitor {
 public:
  SecondPassVisitor(const FirstPassVisitor& first_pass,
                    CPModelProto* const model_proto)
      : expression_map_(first_pass.expression_map()),
        interval_map_(first_pass.interval_map()),
        sequence_map_(first_pass.sequence_map()),
        delegate_map_(first_pass.delegate_map()),
        expression_list_(first_pass.expression_list()),
        constraint_list_(first_pass.constraint_list()),
        interval_list_(first_pass.interval_list()),
        sequence_list_(first_pass.sequence_list()),
        model_proto_(model_proto) {}

  virtual ~SecondPassVisitor() {}

  virtual std::string DebugString() const { return "SecondPassVisitor"; }

  virtual void BeginVisitModel(const std::string& model_name) {
    model_proto_->set_model(model_name);
    model_proto_->set_version(kModelVersion);
    PushArgumentHolder();
    for (const IntExpr* const expr : expression_list_) {
      expr->Accept(this);
    }

    for (const IntervalVar* const var : interval_list_) {
      var->Accept(this);
    }

    for (const SequenceVar* const seq : sequence_list_) {
      seq->Accept(this);
    }
  }

  virtual void EndVisitModel(const std::string& model_name) {
    for (ArgumentHolder* const arg : extensions_) {
      WriteModelExtension(arg);
    }
    PopArgumentHolder();
    // Write tags.
    for (int i = 0; i < tags_.size(); ++i) {
      model_proto_->add_tags(tags_.Element(i));
    }
  }

  virtual void BeginVisitConstraint(const std::string& type_name,
                                    const Constraint* const constraint) {
    PushArgumentHolder();
  }

  virtual void EndVisitConstraint(const std::string& type_name,
                                  const Constraint* const constraint) {
    // We ignore cast constraints, they will be regenerated automatically.
    if (constraint->IsCastConstraint()) {
      return;
    }

    const int index = model_proto_->constraints_size();
    CPConstraintProto* const constraint_proto = model_proto_->add_constraints();
    ExportToProto(constraint, constraint_proto, type_name, index);
    if (constraint->HasName()) {
      constraint_proto->set_name(constraint->name());
    }
    PopArgumentHolder();
  }

  virtual void BeginVisitIntegerExpression(const std::string& type_name,
                                           const IntExpr* const expression) {
    PushArgumentHolder();
  }

  virtual void EndVisitIntegerExpression(const std::string& type_name,
                                         const IntExpr* const expression) {
    const int index = model_proto_->expressions_size();
    CPIntegerExpressionProto* const expression_proto =
        model_proto_->add_expressions();
    ExportToProto(expression, expression_proto, type_name, index);
    PopArgumentHolder();
  }

  virtual void BeginVisitExtension(const std::string& type_name) {
    PushExtension(type_name);
  }

  virtual void EndVisitExtension(const std::string& type_name) {
    PopAndSaveExtension();
  }

  virtual void VisitIntegerArgument(const std::string& arg_name, int64 value) {
    top()->set_integer_argument(arg_name, value);
  }

  virtual void VisitIntegerArrayArgument(const std::string& arg_name,
                                         const std::vector<int64>& values) {
    top()->set_integer_array_argument(arg_name, values);
  }

  virtual void VisitIntegerMatrixArgument(const std::string& arg_name,
                                          const IntTupleSet& values) {
    top()->set_integer_matrix_argument(arg_name, values);
  }

  virtual void VisitIntegerExpressionArgument(const std::string& arg_name,
                                              IntExpr* const argument) {
    top()->set_integer_expression_argument(arg_name,
                                           FindExpressionIndexOrDie(argument));
  }

  virtual void VisitIntegerVariableArrayArgument(
      const std::string& arg_name, const std::vector<IntVar*>& arguments) {
    std::vector<int> indices;
    for (int i = 0; i < arguments.size(); ++i) {
      indices.push_back(FindExpressionIndexOrDie(arguments[i]));
    }
    top()->set_integer_variable_array_argument(arg_name, indices.data(),
                                               indices.size());
  }

  virtual void VisitIntervalArgument(const std::string& arg_name,
                                     IntervalVar* argument) {
    top()->set_interval_argument(arg_name, FindIntervalIndexOrDie(argument));
  }

  virtual void VisitIntervalArrayArgument(
      const std::string& arg_name, const std::vector<IntervalVar*>& arguments) {
    std::vector<int> indices;
    for (int i = 0; i < arguments.size(); ++i) {
      indices.push_back(FindIntervalIndexOrDie(arguments[i]));
    }
    top()->set_interval_array_argument(arg_name, indices.data(),
                                       indices.size());
  }

  virtual void VisitSequenceArgument(const std::string& arg_name,
                                     SequenceVar* argument) {
    top()->set_sequence_argument(arg_name, FindSequenceIndexOrDie(argument));
  }

  virtual void VisitSequenceArrayArgument(
      const std::string& arg_name, const std::vector<SequenceVar*>& arguments) {
    std::vector<int> indices;
    for (int i = 0; i < arguments.size(); ++i) {
      indices.push_back(FindSequenceIndexOrDie(arguments[i]));
    }
    top()->set_sequence_array_argument(arg_name, indices.data(),
                                       indices.size());
  }

  virtual void VisitIntegerVariable(const IntVar* const variable,
                                    IntExpr* const delegate) {
    if (delegate != nullptr) {
      const int index = model_proto_->expressions_size();
      CPIntegerExpressionProto* const var_proto =
          model_proto_->add_expressions();
      var_proto->set_index(index);
      var_proto->set_type_index(TagIndex(ModelVisitor::kIntegerVariable));
      /*
       * MODIFICATION: Add names for cast variables.
       */
      if (variable->HasName()) {
        var_proto->set_name(variable->name());
      }
      /*
       * MODIFICATION: Add names for cast variables.
       */
      CPArgumentProto* const sub_proto = var_proto->add_arguments();
      sub_proto->set_argument_index(
          TagIndex(ModelVisitor::kExpressionArgument));
      sub_proto->set_integer_expression_index(
          FindExpressionIndexOrDie(delegate));
    } else {
      const int index = model_proto_->expressions_size();
      CPIntegerExpressionProto* const var_proto =
          model_proto_->add_expressions();
      var_proto->set_index(index);
      var_proto->set_type_index(TagIndex(ModelVisitor::kIntegerVariable));
      if (variable->HasName()) {
        var_proto->set_name(variable->name());
      }
      if (variable->Size() == variable->Max() - variable->Min() + 1) {
        // Contiguous
        CPArgumentProto* const min_proto = var_proto->add_arguments();
        min_proto->set_argument_index(TagIndex(ModelVisitor::kMinArgument));
        min_proto->set_integer_value(variable->Min());
        CPArgumentProto* const max_proto = var_proto->add_arguments();
        max_proto->set_argument_index(TagIndex(ModelVisitor::kMaxArgument));
        max_proto->set_integer_value(variable->Max());
      } else {
        // Non Contiguous
        CPArgumentProto* const values_proto = var_proto->add_arguments();
        values_proto->set_argument_index(
            TagIndex(ModelVisitor::kValuesArgument));
        std::unique_ptr<IntVarIterator> it(variable->MakeDomainIterator(false));
        for (const int64 value : InitAndGetValues(it.get())) {
          values_proto->add_integer_array(value);
        }
      }
    }
  }

  virtual void VisitIntegerVariable(const IntVar* const variable,
                                    const std::string& operation, int64 value,
                                    IntVar* const delegate) {
    const int index = model_proto_->expressions_size();
    CPIntegerExpressionProto* const var_proto = model_proto_->add_expressions();
    var_proto->set_index(index);
    var_proto->set_type_index(TagIndex(ModelVisitor::kIntegerVariable));
    /*
     * MODIFICATION: Add names for optimized (e.g., x+c, x*c) variables.
     */
    if (variable->HasName()) {
      var_proto->set_name(variable->name());
    }
    /*
     * MODIFICATION: Add names for optimized (e.g., x+c, x*c) variables.
     */
    CPArgumentProto* const sub_proto = var_proto->add_arguments();
    sub_proto->set_argument_index(TagIndex(ModelVisitor::kVariableArgument));
    sub_proto->set_integer_expression_index(FindExpressionIndexOrDie(delegate));
    CPArgumentProto* const value_proto = var_proto->add_arguments();
    value_proto->set_argument_index(TagIndex(operation));
    value_proto->set_integer_value(value);
  }

  virtual void VisitIntervalVariable(const IntervalVar* const variable,
                                     const std::string& operation, int64 value,
                                     IntervalVar* const delegate) {
    if (delegate != nullptr) {
      const int index = model_proto_->intervals_size();
      CPIntervalVariableProto* const var_proto = model_proto_->add_intervals();
      var_proto->set_index(index);
      var_proto->set_type_index(TagIndex(ModelVisitor::kIntervalVariable));
      CPArgumentProto* const sub_proto = var_proto->add_arguments();
      sub_proto->set_argument_index(TagIndex(operation));
      sub_proto->set_interval_index(FindIntervalIndexOrDie(delegate));
      sub_proto->set_integer_value(value);
      if (operation == ModelVisitor::kStartSyncOnStartOperation ||
          operation == ModelVisitor::kStartSyncOnEndOperation) {
        CHECK_EQ(delegate->DurationMin(), delegate->DurationMax());
        sub_proto->add_integer_array(delegate->DurationMin());
      }
    } else {
      const int index = model_proto_->intervals_size();
      CPIntervalVariableProto* const var_proto = model_proto_->add_intervals();
      var_proto->set_index(index);
      var_proto->set_type_index(TagIndex(ModelVisitor::kIntervalVariable));
      if (variable->HasName()) {
        var_proto->set_name(variable->name());
      }
      CPArgumentProto* const start_min_proto = var_proto->add_arguments();
      start_min_proto->set_argument_index(
          TagIndex(ModelVisitor::kStartMinArgument));
      start_min_proto->set_integer_value(variable->StartMin());
      CPArgumentProto* const start_max_proto = var_proto->add_arguments();
      start_max_proto->set_argument_index(
          TagIndex(ModelVisitor::kStartMaxArgument));
      start_max_proto->set_integer_value(variable->StartMax());
      CPArgumentProto* const end_min_proto = var_proto->add_arguments();
      end_min_proto->set_argument_index(
          TagIndex(ModelVisitor::kEndMinArgument));
      end_min_proto->set_integer_value(variable->EndMin());
      CPArgumentProto* const end_max_proto = var_proto->add_arguments();
      end_max_proto->set_argument_index(
          TagIndex(ModelVisitor::kEndMaxArgument));
      end_max_proto->set_integer_value(variable->EndMax());
      CPArgumentProto* const duration_min_proto = var_proto->add_arguments();
      duration_min_proto->set_argument_index(
          TagIndex(ModelVisitor::kDurationMinArgument));
      duration_min_proto->set_integer_value(variable->DurationMin());
      CPArgumentProto* const duration_max_proto = var_proto->add_arguments();
      duration_max_proto->set_argument_index(
          TagIndex(ModelVisitor::kDurationMaxArgument));
      duration_max_proto->set_integer_value(variable->DurationMax());
      CPArgumentProto* const optional_proto = var_proto->add_arguments();
      optional_proto->set_argument_index(
          TagIndex(ModelVisitor::kOptionalArgument));
      optional_proto->set_integer_value(!variable->MustBePerformed());
    }
  }

  virtual void VisitSequenceVariable(const SequenceVar* const sequence) {
    const int index = model_proto_->sequences_size();
    CPSequenceVariableProto* const var_proto = model_proto_->add_sequences();
    var_proto->set_index(index);
    var_proto->set_type_index(TagIndex(ModelVisitor::kSequenceVariable));
    if (sequence->HasName()) {
      var_proto->set_name(sequence->name());
    }
    CPArgumentProto* const sub_proto = var_proto->add_arguments();
    sub_proto->set_argument_index(TagIndex(ModelVisitor::kIntervalsArgument));
    for (int i = 0; i < sequence->size(); ++i) {
      IntervalVar* const interval = sequence->Interval(i);
      sub_proto->add_interval_array(FindIntervalIndexOrDie(interval));
    }
  }

  int TagIndex(const std::string& tag) { return tags_.Add(tag); }

 private:
  void WriteModelExtension(ArgumentHolder* const holder) {
    CHECK(holder != nullptr);
    if (holder->type_name().compare(kObjectiveExtension) == 0) {
      WriteObjective(holder);
    } else if (holder->type_name().compare(kSearchLimitExtension) == 0) {
      WriteSearchLimit(holder);
    } else if (holder->type_name().compare(kVariableGroupExtension) == 0) {
      WriteVariableGroup(holder);
    } else {
      LOG(INFO) << "Unknown model extension :" << holder->type_name();
    }
  }

  void WriteObjective(ArgumentHolder* const holder) {
    CHECK(holder != nullptr);
    const bool maximize = holder->FindIntegerArgumentOrDie(kMaximizeArgument);
    const int64 step = holder->FindIntegerArgumentOrDie(kStepArgument);
    const int objective_index =
        holder->FindIntegerExpressionArgumentOrDie(kExpressionArgument);
    CPObjectiveProto* const objective_proto = model_proto_->mutable_objective();
    objective_proto->set_maximize(maximize);
    objective_proto->set_step(step);
    objective_proto->set_objective_index(objective_index);
  }

  void WriteSearchLimit(ArgumentHolder* const holder) {
    CHECK(holder != nullptr);
    SearchLimitProto* const proto = model_proto_->mutable_search_limit();
    proto->set_time(
        holder->FindIntegerArgumentWithDefault(kTimeLimitArgument, kint64max));
    proto->set_branches(holder->FindIntegerArgumentWithDefault(
        kBranchesLimitArgument, kint64max));
    proto->set_failures(holder->FindIntegerArgumentWithDefault(
        kFailuresLimitArgument, kint64max));
    proto->set_solutions(holder->FindIntegerArgumentWithDefault(
        kSolutionLimitArgument, kint64max));
    proto->set_smart_time_check(
        holder->FindIntegerArgumentWithDefault(kSmartTimeCheckArgument, false));
    proto->set_cumulative(
        holder->FindIntegerArgumentWithDefault(kCumulativeArgument, false));
  }

  void WriteVariableGroup(ArgumentHolder* const holder) {
    CPVariableGroup* const group_proto = model_proto_->add_variable_groups();
    holder->ExportToProto(&tags_, group_proto);
  }

  template <class A, class P>
  void ExportToProto(const A* const argument, P* const proto,
                     const std::string& type_name, int index) {
    CHECK(proto != nullptr);
    CHECK(argument != nullptr);
    proto->set_index(index);
    proto->set_type_index(TagIndex(type_name));
    if (argument->HasName()) {
      proto->set_name(argument->name());
    }
    top()->ExportToProto(&tags_, proto);
    for (ArgumentHolder* const arg : extensions_) {
      CPExtensionProto* const extension_proto = proto->add_extensions();
      extension_proto->set_type_index(TagIndex(arg->type_name()));
      arg->ExportToProto(&tags_, extension_proto);
    }
  }

  void PushArgumentHolder() { holders_.push_back(new ArgumentHolder); }

  void PopArgumentHolder() {
    CHECK(!holders_.empty());
    delete holders_.back();
    holders_.pop_back();
    STLDeleteElements(&extensions_);
    extensions_.clear();
  }

  void PushExtension(const std::string& type_name) {
    PushArgumentHolder();
    holders_.back()->set_type_name(type_name);
  }

  void PopAndSaveExtension() {
    CHECK(!holders_.empty());
    extensions_.push_back(holders_.back());
    holders_.pop_back();
  }

  ArgumentHolder* top() const {
    CHECK(!holders_.empty());
    return holders_.back();
  }

  int FindExpressionIndexOrDie(IntExpr* const expression) const {
    return FindOrDie(expression_map_, expression);
  }

  int FindIntervalIndexOrDie(IntervalVar* const interval) const {
    return FindOrDie(interval_map_, interval);
  }

  int FindSequenceIndexOrDie(SequenceVar* const sequence) const {
    return FindOrDie(sequence_map_, sequence);
  }

  hash_map<const IntExpr*, int> expression_map_;
  hash_map<const IntervalVar*, int> interval_map_;
  hash_map<const SequenceVar*, int> sequence_map_;
  hash_map<const IntVar*, IntExpr*> delegate_map_;
  std::vector<const IntExpr*> expression_list_;
  std::vector<const Constraint*> constraint_list_;
  std::vector<const IntervalVar*> interval_list_;
  std::vector<const SequenceVar*> sequence_list_;
  CPModelProto* const model_proto_;

  std::vector<ArgumentHolder*> holders_;
  std::vector<ArgumentHolder*> extensions_;
  VectorMap<std::string> tags_;
};
} /* namespace <anonymous> */
} /* namespace operations_research */

/*
 * This is the or-tools loader class for the model. It loads the model from
 * a protobuf, to which the model was previously exported by visitors.
 *
 * We need it to properly load the model. We register our own builder for
 * UDF functions and use some of the loader's methods. Thus, we need at least
 * partial definition. The rest will be resolved dynamically by the linker.
 *
 * Since the template definition is given here, the function will be defined
 * HERE, not in the constraint solver library. This means any changes in
 * the definition of the CPModelLoader class in the library might cause
 * binary compatibility issues and nasty bugs. Should be fine for this project
 * though.
 */
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

namespace searchlight {

namespace {

using operations_research::CPIntegerExpressionProto;
using operations_research::CPModelLoader;
using operations_research::FirstPassVisitor;
using operations_research::SecondPassVisitor;
using operations_research::CPConstraintProto;

/*
 * Builds UDFs during model duplication (at import).
 */
IntExpr* UDFBuilder(UDFFunctionCreator udf_creator, const AdapterPtr &adapter,
        CPModelLoader * const builder,
        const CPIntegerExpressionProto &proto) {

    // should have IntVar parameters protobuffed
    std::vector<IntVar *> vars;
    if (!builder->ScanArguments(ModelVisitor::kVarsArgument, proto, &vars)) {
        return nullptr;
    }

    // should have integer parameters protobuffed
    std::vector<int64> params;
    if (!builder->ScanArguments(ModelVisitor::kValuesArgument, proto,
            &params)) {
        return nullptr;
    }

    return builder->solver()->RevAlloc(
            udf_creator(builder->solver(), adapter, vars, params));
}

// Relaxable constraint builders
Constraint *BuildRelaxableLessOrEqual(CPModelLoader* const builder,
                             const CPConstraintProto& proto) {
    int64 id;
    if (!builder->ScanArguments(RelaxableConstraint::ModelIDTag, proto, &id)) {
    	return nullptr;
    }
    IntExpr *expr;
    if (!builder->ScanArguments(ModelVisitor::kExpressionArgument, proto,
    		&expr)) {
    	return nullptr;
    }
	int64 value;
	if (!builder->ScanArguments(ModelVisitor::kValueArgument, proto,
			&value)) {
		return nullptr;
	}
	// All is valid; create
	Solver *solver = builder->solver();
	RelaxableConstraint *constr = solver->RevAlloc(
			new LessEqExprCst(solver, expr, value));
	constr->SetId(id);
	return constr;
}

Constraint *BuildRelaxableGreaterOrEqual(CPModelLoader* const builder,
                             const CPConstraintProto& proto) {
    int64 id;
    if (!builder->ScanArguments(RelaxableConstraint::ModelIDTag, proto, &id)) {
    	return nullptr;
    }
    IntExpr *expr;
    if (!builder->ScanArguments(ModelVisitor::kExpressionArgument, proto,
    		&expr)) {
    	return nullptr;
    }
	int64 value;
	if (!builder->ScanArguments(ModelVisitor::kValueArgument, proto,
			&value)) {
		return nullptr;
	}
	// All is valid; create
	Solver *solver = builder->solver();
	RelaxableConstraint *constr = solver->RevAlloc(
			new GreaterEqExprCst(solver, expr, value));
	constr->SetId(id);
	return constr;
}

Constraint *BuildRelaxableBetween(CPModelLoader* const builder,
                             const CPConstraintProto& proto) {
    int64 id;
    if (!builder->ScanArguments(RelaxableConstraint::ModelIDTag, proto, &id)) {
    	return nullptr;
    }
	int64 min_value;
	if (!builder->ScanArguments(ModelVisitor::kMinArgument, proto,
			&min_value)) {
		return nullptr;
	}
    IntExpr *expr;
    if (!builder->ScanArguments(ModelVisitor::kExpressionArgument, proto,
    		&expr)) {
    	return nullptr;
    }
	int64 max_value;
	if (!builder->ScanArguments(ModelVisitor::kMaxArgument, proto,
			&max_value)) {
		return nullptr;
	}
	// All is valid; create
	Solver *solver = builder->solver();
	RelaxableConstraint *constr = solver->RevAlloc(
			new BetweenCt(solver, expr, min_value, max_value));
	constr->SetId(id);
	return constr;
}

// Register builders for relaxable constraints
bool RegisterRelaxableCstBuilders(Solver &solver) {
	solver.RegisterBuilder(RelaxableConstraint::BetweenConstTag,
			NewPermanentCallback(BuildRelaxableBetween));
	solver.RegisterBuilder(RelaxableConstraint::LessEqConstTag,
			NewPermanentCallback(BuildRelaxableLessOrEqual));
	solver.RegisterBuilder(RelaxableConstraint::GreaterEqConstTag,
			NewPermanentCallback(BuildRelaxableGreaterOrEqual));
	return true;
}

/*
 * Registers a builder for UDFs with the destination solver.
 */
bool RegisterUDFBuilder(const Searchlight &sl, Solver &solver,
        const AdapterPtr &adapter) {
    StringSet udfs = sl.GetAllUsedUDFs();
    for (StringSet::const_iterator cit = udfs.begin(); cit != udfs.end();
            cit++) {
        UDFFunctionCreator udf_creator = sl.GetRegisteredUDF(*cit);
        if (!udf_creator) {
            return false;
        }
        solver.RegisterBuilder(*cit,
                NewPermanentCallback(UDFBuilder, udf_creator, adapter));
        /*
         * The callback (an object) is deleted by the solver at dtor.
         */
    }
    return true;
}
} /* namespace <anonymous> */

void ExportModel(const Solver &solver, CPModelProto *model_proto) {
  FirstPassVisitor first_pass;
  solver.Accept(&first_pass);
  SecondPassVisitor second_pass(first_pass, model_proto);
  solver.Accept(&second_pass);
}

bool CloneModel(const Searchlight &sl, const Solver &from, Solver &to,
        const AdapterPtr &adapter) {
    // Export the model
    CPModelProto model;
    ExportModel(from, &model);

    // Clone
    return CloneModel(sl, model, to, adapter);
}

bool CloneModel(const Searchlight &sl, const CPModelProto &from, Solver &to,
        const AdapterPtr &adapter) {

    // Register UDF and relaxable constraints builder
    if (!RegisterUDFBuilder(sl, to, adapter) ||
    		!RegisterRelaxableCstBuilders(to)) {
        return false;
    }

    // Import the model
    return to.LoadModel(from);
}

} /* namespace searchlight */
