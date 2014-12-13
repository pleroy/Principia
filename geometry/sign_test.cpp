#include "geometry/sign.hpp"

#include "base/heap_checker.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace principia {
namespace geometry {

class SignTest : public testing::Test {
 protected:
  Sign const positive_ = Sign(1);
  Sign const negative_ = Sign(-1);

 private:
  base::HeapChecker heap_checker_;
};

TEST_F(SignTest, Integer) {
  EXPECT_TRUE(positive_.Positive());
  EXPECT_FALSE(positive_.Negative());
  EXPECT_FALSE(negative_.Positive());
  EXPECT_TRUE(negative_.Negative());
}

TEST_F(SignTest, SignMultiplication) {
  EXPECT_TRUE((positive_ * positive_).Positive());
  EXPECT_TRUE((positive_ * negative_).Negative());
  EXPECT_TRUE((negative_ * positive_).Negative());
  EXPECT_TRUE((negative_ * negative_).Positive());
}

TEST_F(SignTest, ScalarMultiplication) {
  EXPECT_EQ(3, positive_ * 3);
  EXPECT_EQ(-3, positive_ * -3);
  EXPECT_EQ(-3, negative_ * 3);
  EXPECT_EQ(3, negative_ * -3);
}

}  // namespace geometry
}  // namespace principia
