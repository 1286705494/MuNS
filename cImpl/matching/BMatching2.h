/*
 * BMatching2.h
 *
 *  Created on: 02/09/2014
 *      Author: jefcha
 */

#ifndef BMATCHING2_H_
#define BMATCHING2_H_


#include <vector>

extern
//float matching(int n, int m, const std::vector<double>& nzv, const std::vector<int>& nzi, const std::vector<int>& nzj, int* m1, int* m2);
double matching2(int n, int m, const std::vector<float>& vA, std::vector<int>& m1, std::vector<int>& m2);





#endif /* BMATCHING2_H_ */
