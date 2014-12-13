﻿#include "testing_utilities/vanishes_before.hpp"

#include <limits>

#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/bipm.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/heap_checked_test.hpp"

using principia::bipm::Knot;
using principia::quantities::Speed;
using testing::Ne;

namespace principia {
namespace testing_utilities {

class VanishesBeforeTest : public testing_utilities::HeapCheckedTest {};

TEST_F(VanishesBeforeTest, Dimensionless) {
  double const y = 3000.0 * std::numeric_limits<double>::epsilon();
  EXPECT_THAT(y, VanishesBefore(1000.0, 6));
  EXPECT_THAT(2 * y, Not(VanishesBefore(1000.0, 6)));
  double const δy = e / 100.0;
  double e_accumulated = 0.0;
  for (int i = 1; i <= 100.0; ++i) {
    e_accumulated += δy;
  }
  EXPECT_THAT(e_accumulated, Ne(e));
  EXPECT_THAT(e_accumulated - e, Not(VanishesBefore(e, 4)));
  EXPECT_THAT(e_accumulated - e, VanishesBefore(e, 1));
}

TEST_F(VanishesBeforeTest, Quantity) {
  Speed v1 = 1 * Knot;
  Speed const v2 = 3 * v1 * std::numeric_limits<double>::epsilon();
  EXPECT_THAT(v2, VanishesBefore(v1, 3));
  EXPECT_THAT(2 * v2, Not(VanishesBefore(v1, 3)));
  Speed const δv = v1 / 100;
  Speed v_accumulated;
  for (int i = 1; i <= 100; ++i) {
    v_accumulated += δv;
  }
  EXPECT_THAT(v_accumulated, Ne(v1));
  EXPECT_THAT(v_accumulated - v1, Not(VanishesBefore(v1, 8)));
  EXPECT_THAT(v_accumulated - v1, VanishesBefore(v1, 4));
}

}  // namespace testing_utilities
}  // namespace principia
