#ifndef LTEST_HPP
#define LTEST_HPP

#include <iostream>
#include <cmath>

#define TEST_EPSILON         1e-7

#define TEST_START(name) std::cout << "TEST " << (name) << " -----------------------------------" << std::endl;
#define TEST_END(name) std::cout << std::endl;

#define EXPECT_NEAR(v1, v2, esp) {\
  if ( std::abs((v1) - (v2)) <= (esp) )\
    std::cout << "EXPECT_NEAR( " << v1 << ", " << v2 << " ) ... " << "[success]" << std::endl;\
  else\
    std::cout << "EXPECT_NEAR( " << v1 << ", " << v2 << " ) ... " << "[failed]" << std::endl;\
}

#define EXPECT_EQ(v1, v2) {\
  if ( (v1) == (v2) )\
    std::cout << "EXPECT_EQ( " << v1 << ", " << v2 << " ) ... " << "[success]" << std::endl;\
  else\
    std::cout << "EXPECT_EQ( " << v1 << ", " << v2 << " ) ... " << "[failed]" << std::endl;\
}

#define EXPECT_LE(v1, v2) {\
  if ( (v1) <= (v2) )\
    std::cout << "EXPECT_LE( " << v1 << ", " << v2 << " ) ... " << "[success]" << std::endl;\
  else\
    std::cout << "EXPECT_LE( " << v1 << ", " << v2 << " ) ... " << "[failed]" << std::endl;\
}

#define EXPECT_GE(v1, v2) {\
  if ( (v1) <= (v2) )\
    std::cout << "EXPECT_GE( " << v1 << ", " << v2 << " ) ... " << "[success]" << std::endl;\
  else\
    std::cout << "EXPECT_GE( " << v1 << ", " << v2 << " ) ... " << "[failed]" << std::endl;\
}

#define EXPECT_DOUBLE_EQ(v1, v2) {\
  if ( std::abs((v1) - (v2)) <= TEST_EPSILON )\
    std::cout << "EXPECT_DOUBLE_EQ( " << v1 << ", " << v2 << " ) ... " << "[success]" << std::endl;\
  else\
    std::cout << "EXPECT_DOUBLE_EQ( " << v1 << ", " << v2 << " ) ... " << "[failed]" << std::endl;\
}

#define EXPECT_FLOAT_EQ(v1, v2) {\
  if ( std::abs((v1) - (v2)) <= TEST_EPSILON )\
    std::cout << "EXPECT_FLOAT_EQ( " << v1 << ", " << v2 << " ) ... " << "[success]" << std::endl;\
  else\
    std::cout << "EXPECT_FLOAT_EQ( " << v1 << ", " << v2 << " ) ... " << "[failed]" << std::endl;\
}

#define EXPECT_FALSE(v1) {\
  if ( (v1) == false )\
    std::cout << "EXPECT_FALSE( " << v1 << " ) ... " << "[success]" << std::endl;\
  else\
    std::cout << "EXPECT_FALSE( " << v1 << " ) ... " << "[failed]" << std::endl;\
}

#define EXPECT_TRUE(v1) {\
  if ( (v1) == true )\
    std::cout << "EXPECT_TRUE( " << v1 << " ) ... " << "[success]" << std::endl;\
  else\
    std::cout << "EXPECT_TRUE( " << v1 << " ) ... " << "[failed]" << std::endl;\
}

#endif