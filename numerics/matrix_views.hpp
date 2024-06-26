#pragma once

#include <concepts>
#include <utility>

#include "numerics/concepts.hpp"
#include "numerics/transposed_view.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _matrix_views {
namespace internal {

using namespace principia::numerics::_concepts;
using namespace principia::numerics::_transposed_view;
using namespace principia::quantities::_named_quantities;

// A view of a rectangular block of a matrix.  This view is |two_dimensional|.
template<typename Matrix>
  requires two_dimensional<Matrix>
struct BlockView {
  using Scalar = typename Matrix::Scalar;

  Matrix& matrix;
  int first_row;
  int last_row;
  int first_column;
  int last_column;

  constexpr Scalar& operator()(int row, int column);
  constexpr Scalar const& operator()(int row, int column) const;

  template<typename T>
    requires two_dimensional<T> && same_elements_as<T, Matrix>
  BlockView& operator=(T const& right);

  template<typename T>
    requires two_dimensional<T> && same_elements_as<T, Matrix>
  BlockView& operator+=(T const& right);
  template<typename T>
    requires two_dimensional<T> && same_elements_as<T, Matrix>
  BlockView& operator-=(T const& right);

  BlockView& operator*=(double right);
  BlockView& operator/=(double right);

  constexpr int rows() const;
  constexpr int columns() const;
};

// A view of a column of a matrix.  This view is |one_dimensional|.
template<typename Matrix>
  requires two_dimensional<Matrix>
struct ColumnView {
  using Scalar = typename Matrix::Scalar;

  Matrix& matrix;
  int first_row;
  int last_row;
  int column;

  constexpr Scalar& operator[](int index);
  constexpr Scalar const& operator[](int index) const;

  template<typename T>
    requires one_dimensional<T> && same_elements_as<T, Matrix>
  ColumnView& operator=(T const& right);
  ColumnView& operator=(ColumnView const& right);

  template<typename T>
    requires one_dimensional<T> && same_elements_as<T, Matrix>
  ColumnView& operator+=(T const& right);
  template<typename T>
    requires one_dimensional<T> && same_elements_as<T, Matrix>
  ColumnView& operator-=(T const& right);

  ColumnView& operator*=(double right);
  ColumnView& operator/=(double right);

  Scalar Norm() const;
  Square<Scalar> Norm²() const;
  constexpr int size() const;
};

// TODO(phl): This should probably be just |swap|.  The semantics of BlockView
// and ColumnView do imply an implicit dereferencing, so swapping should work
// the same.
template<typename Matrix>
void SwapColumns(ColumnView<Matrix>& m1, ColumnView<Matrix>& m2);

template<typename LMatrix, typename RMatrix>
Product<typename LMatrix::Scalar, typename RMatrix::Scalar> operator*(
    TransposedView<ColumnView<LMatrix>> const& left,
    ColumnView<RMatrix> const& right);

// The declarations below also declare the symmetric operator== and
// operator!=.

template<typename Matrix, typename T>
  requires two_dimensional<T> && same_elements_as<T, Matrix>
bool operator==(BlockView<Matrix> const& left, T const& right);

template<typename Matrix, typename T>
  requires one_dimensional<T> && same_elements_as<T, Matrix>
bool operator==(ColumnView<Matrix> const& left,
                T const& right);

template<typename Matrix>
std::ostream& operator<<(std::ostream& out,
                         BlockView<Matrix> const& view);

template<typename Matrix>
std::ostream& operator<<(std::ostream& out,
                         ColumnView<Matrix> const& view);

}  // namespace internal

using internal::BlockView;
using internal::ColumnView;
using internal::SwapColumns;

}  // namespace _matrix_views
}  // namespace numerics
}  // namespace principia

#include "numerics/matrix_views_body.hpp"
