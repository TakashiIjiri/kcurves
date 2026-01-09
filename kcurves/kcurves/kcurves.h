#pragma once

#pragma unmanaged 

//Eigen 
#include "./3rdparty/Eigen/Dense"
#include "./3rdparty/Eigen/Geometry"
#include "./3rdparty/Eigen/Sparse"

typedef Eigen::Vector2i EVec2i;
typedef Eigen::Vector2f EVec2f;
typedef Eigen::SparseMatrix<float> ESpMat;
typedef Eigen::Triplet<float>      ETrip;

#include <vector>
#include <iostream>


inline void Trace(EVec2f& p)
{
	std::cout << p[0] << "" << p[1] << std::endl;
}




//実数解を一つ返す
//solve x^3 + a x^2 + b x + c
inline float SolveCubicFunction
(
	const float a,
	const float b,
	const float c
)
{
	const float p = b - a * a /  3.0f;
	const float q = 2 * a*a*a / 27.0f - a * b / 3.0f + c;

	const float tmp = q * q / 4.0f + p * p * p / 27.0f;
	if( tmp < 0) return 0;

	const float mTmp = - q / 2.0f + sqrt(tmp);
	const float nTmp = - q / 2.0f - sqrt(tmp);
	const float m = (mTmp < 0) ?  - pow( - mTmp, 1/3.0f) :  pow(mTmp, 1/3.0f);
	const float n = (nTmp < 0) ?  - pow( - nTmp, 1/3.0f) :  pow(nTmp, 1/3.0f);

	return - a / 3.0f + m + n;
}





inline float TriangleArea(
		const EVec2f &p0, 
		const EVec2f &p1, 
		const EVec2f &p2)
{
	EVec2f v1 = p1-p0;
	EVec2f v2 = p2-p0;
	const float cross = v1[0] * v2[1] - v1[1] * v2[0];
	return 0.5f * std::abs(cross);
}





// input : cps (2D control points,        size:N)
// output: out_cps  (2D Bezier control points, size:3 * N, each 3 construct quad bezier)
// output:λi  (weight to compute Bezier cp, size:N)
//
// note
// 各制御点 cps[i] に対して、2次ベジエ制御点 Ci0, Ci1, Ci2を生成する
//
// ただし，
// C_i1 は上記の C1[i]に対応
// C_i0, C_i+1,0 については、『C_i0 = C_i+1,0 = (1-λi) Ci + λi Ci+1 』と計算できるのでλiとC1[i]のみ保持する
//
//https://qiita.com/Rijicho_nl/items/05ee4c8d77e99e29daa5
//
const int ITER_NUM = 15;	

inline void compute_kCurves
(
	const std::vector<EVec2f> &cps,

	int  num_sample, //sampling number of each curve segment
	std::vector<EVec2f> &out_cps,
	std::vector<EVec2f> &out_points
)
{
	if( cps.size() < 4 ) return;
	if( num_sample < 3) num_sample = 3;

	const int N = (int)cps.size();
	std::vector<float> lambda(N, 0.5f), ti(N);
	std::vector<EVec2f> Ci0(N), Ci1 = cps, Ci2(N);

	for( int iter = 0; iter < ITER_NUM; ++iter)
	{
		//1. update lambda (use 0.5 for the first iteration (更新してもいい気もするけど..))
		if( iter != 0 ) 
		{
			for( int i = 0; i < N; ++i)
			{
				int next_i = (i+1) % N;
				float A = sqrt( TriangleArea(Ci0[i], Ci1[i     ], Ci1[next_i] ) );
				float B = sqrt( TriangleArea(Ci0[i], Ci1[i     ], Ci1[next_i] ) );
				float C = sqrt( TriangleArea(Ci1[i], Ci1[next_i], Ci2[next_i] ) );
				lambda[i] = A / (B + C);
			}
		}

		//2. update Ci2, Ci+10
		for( int i = 0; i < N; ++i)
		{
			int next_i = (i+1)%N;
			Ci2[i] = Ci0[next_i] = (1 - lambda[i]) * Ci1[  i   ]  +  
				                          lambda[i]  * Ci1[next_i];
		}

		//3. update ti
		for( int i = 0; i < N; ++i)
		{
			//solve eq(4)
			EVec2f Ci2_Ci0 = Ci2[i] - Ci0[i];
			EVec2f Ci0_pi  = Ci0[i] - cps[i];
			float a = Ci2_Ci0.squaredNorm();
			float b = 3 * Ci2_Ci0.dot( Ci0_pi );
			float c = (3 * Ci0[i] - 2 * cps[i] - Ci2[i]).dot( Ci0_pi );
			float d = - Ci0_pi.squaredNorm();

			ti[i] = (a==0) ? 0.5f : SolveCubicFunction(b/a, c/a, d/a);
		}

		//4. update C1
		ESpMat A(N,N);
		Eigen::VectorXf b1(N), b2(N);
		std::vector< Eigen::Triplet<float> > entries; //{row, col, val}

		for( int i = 0; i < N; ++i)
		{
			int nextI = (i+1) % N;
			int prevI = (i-1) < 0 ? N - 1 : i - 1;
			float t = ti[i];
			float a = (1-lambda[prevI]) * (1 - t) * (1 - t);
			float b = lambda[i] * t * t;
			float c = lambda[prevI] * (1 - t) * (1 - t) + (2-(1+lambda[i]) * t) * t;

			entries.push_back( ETrip(i, prevI, a) );
			entries.push_back( ETrip(i, i    , c) );
			entries.push_back( ETrip(i, nextI, b) );
			b1[i] = cps[i][0];
			b2[i] = cps[i][1];
		}

		A.setFromTriplets(entries.begin(), entries.end());
		// solve Ax = b
		Eigen::SparseLU<ESpMat> LU(A);  
		Eigen::VectorXf x = LU.solve(b1);
		Eigen::VectorXf y = LU.solve(b2);

		for( int i = 0; i < N; ++i) Ci1[i] << x[i], y[i];

	}

	out_cps.clear();
	for( int i = 0; i < N; ++i)
	{
		out_cps.push_back( Ci0[i] );
		out_cps.push_back( Ci1[i] );
		out_cps.push_back( Ci2[i] );
	}

	out_points.clear();
	for( int i = 0; i < N; ++i)
	{
		const EVec2f &C0 = Ci0[i];
		const EVec2f &C1 = Ci1[i];
		const EVec2f &C2 = Ci2[i];

		for( int j=0; j < num_sample; ++j)
		{
			float t = j * 1.0f / num_sample;
			EVec2f p = (1 - t) * (1 - t) * C0 + 2 * t * (1 - t) * C1 + t * t * C2;
			out_points.push_back( p );
		}
	}
}





/*
//   solve 
//   2  3  0  0  0     x1       8 
//   3  0  4  0  6     x2      45
//   0 -1 -3  2  0     x3    = -3
//   0  0  1  0  0     x4       3
//   0  4  2  0  1     x5      19
inline void EigenSparseMatPractice()
{
	//prepare field 
	ESpMat A(5, 5);
	Eigen::VectorXd b(5);

	//fill A
	vector< Eigen::Triplet<double> > entries; //{row, col, val}
	entries.push_back(Eigen::Triplet<double>(0, 0, 2));
	entries.push_back(Eigen::Triplet<double>(1, 0, 3));

	entries.push_back(Eigen::Triplet<double>(0, 1, 3));
	entries.push_back(Eigen::Triplet<double>(2, 1, -1));
	entries.push_back(Eigen::Triplet<double>(4, 1, 4));

	entries.push_back(Eigen::Triplet<double>(2, 3, 2));

	entries.push_back(Eigen::Triplet<double>(1, 4, 6));
	entries.push_back(Eigen::Triplet<double>(4, 4, 1));

	entries.push_back(Eigen::Triplet<double>(1, 2, 4));
	entries.push_back(Eigen::Triplet<double>(2, 2, -3));
	entries.push_back(Eigen::Triplet<double>(3, 2, 1));
	entries.push_back(Eigen::Triplet<double>(4, 2, 2));


	A.setFromTriplets(entries.begin(), entries.end());

	// fill b
	b[0] = 8;
	b[1] = 45;
	b[2] = -3;
	b[3] = 3;
	b[4] = 19;

	// solve Ax = b
	Eigen::SparseLU<ESpMat> LU(A);
	Eigen::VectorXd x = LU.solve(b);

	printf("%f %f %f %f %f\n", x[0], x[1], x[2], x[3], x[4]);
	return;
}
*/





#pragma managed 



