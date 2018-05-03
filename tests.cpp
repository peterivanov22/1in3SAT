#include "solver.cpp"
#include <math.h>

#include <glpk.h>
#include "solver.h"

#include <time.h>
#include <stdlib.h>

#include <vector>
#include <list>
#include <utility>
#include <cstdio>

#include <queue>
#include <iostream>
#include <set>
#include <algorithm>

#include <cmath>
#include <gtest/gtest.h>

int argc; // Making arg and arv global to access within TESTs
char ** argv;

TEST(SolverTest, correct1) { 
	int solutions[5] = {0,1,0,0,1};
	int * output = solve(argc,argv, 1);

	ASSERT_EQ(solutions[0], output[0]);
	ASSERT_EQ(solutions[1], output[1]);
	ASSERT_EQ(solutions[2], output[2]);
	ASSERT_EQ(solutions[3], output[3]);
	ASSERT_EQ(solutions[4], output[4]);




}

TEST(SolverTest, correct2) { 
	int * output2 = solve(argc,argv, 2);

	ASSERT_EQ(0, output2[0]);
	ASSERT_EQ(0, output2[1]);
	ASSERT_EQ(0, output2[2]);
	ASSERT_EQ(1, output2[3]);
	ASSERT_EQ(1, output2[4]);


}



int main(int t_argc, char **t_argv) {

    argc = t_argc;
    argv = t_argv;
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}