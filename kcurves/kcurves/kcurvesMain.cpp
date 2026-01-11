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
  MainForm::getInst()->ShowDialog();
  return 0;
}
