// kcurves.cpp : メイン プロジェクト ファイルです。

#include "stdafx.h"
#include "MainForm.h"

#include "kcurves.h"

using namespace System;
using namespace kcurves;

#pragma managed

int main()
{
    Console::WriteLine(L"Hello World");


	EigenSparseMatPractice();

	double x = getRealSolutionOfCubicFunc(1,2,3);
	printf( "solved %f %f\n", x, x*x*x + 1*x*x + 2*x + 3);
	x = getRealSolutionOfCubicFunc(3,2,1);
	printf( "solved %f %f\n", x, x*x*x + 3*x*x + 2*x + 1);
	x = getRealSolutionOfCubicFunc(2,5,1);
	printf( "solved %f %f\n", x, x*x*x + 2*x*x + 5*x + 1);
	x = getRealSolutionOfCubicFunc(1,2,-4);
	printf( "solved %f %f\n", x, x*x*x + 1*x*x + 2*x - 4);




	MainForm::getInst()->ShowDialog();

    return 0;
}
