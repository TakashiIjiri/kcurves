#include "stdafx.h"
#include "MainForm.h"

#include "KCurveUI.h"


using namespace kcurves;


System::Void MainForm::MainForm_MouseUp(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)  
{
	if (e->Button == System::Windows::Forms::MouseButtons::Left  ) KCurveUI::getInst()->LBtnUp( e->X, e->Y);
	if (e->Button == System::Windows::Forms::MouseButtons::Middle) KCurveUI::getInst()->MBtnUp( e->X, e->Y);
	if (e->Button == System::Windows::Forms::MouseButtons::Right ) KCurveUI::getInst()->RBtnUp( e->X, e->Y);
}

System::Void MainForm::MainForm_MouseMove(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
	KCurveUI::getInst()->MouseMove( e->X, e->Y);
}

System::Void MainForm::MainForm_MouseDown(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
	if (e->Button == System::Windows::Forms::MouseButtons::Left  ) KCurveUI::getInst()->LBtnDown( e->X, e->Y);
	if (e->Button == System::Windows::Forms::MouseButtons::Middle) KCurveUI::getInst()->MBtnDown( e->X, e->Y);
	if (e->Button == System::Windows::Forms::MouseButtons::Right ) KCurveUI::getInst()->RBtnDown( e->X, e->Y);
}






void MainForm::RepaintFunction(Object^ sender, PaintEventArgs^ e)
{
	const vector<EVec2d> &CPs      = KCurveUI::getInst()->m_CPs;
	const vector<EVec2d> &kCurveCP = KCurveUI::getInst()->m_kCurveCP;
	const vector<EVec2d> &points   = KCurveUI::getInst()->m_curves;

	System::Drawing::Graphics^  g = e->Graphics;

	for( const auto &p : CPs)  g->DrawEllipse(gcnew Pen(Color::Red,3), (int)p[0]-CIRCLE_R, (int)p[1]-CIRCLE_R, CIRCLE_R*2, CIRCLE_R*2);

	
	for( int i=0; i < (int)kCurveCP.size(); ++i)
	{
		int nexI = (i+1) % (int)kCurveCP.size();
		g->DrawLine(gcnew Pen(Color::LightBlue,2), 
			(int)kCurveCP[i][0], (int)kCurveCP[i][1], 
			(int)kCurveCP[nexI][0], (int)kCurveCP[nexI][1] );
	}


	for( int i=0; i < (int)points.size(); ++i)
	{
		int nexI = (i+1) % (int)points.size();
		g->DrawLine(gcnew Pen(Color::Blue,3), 
			(int)points[i   ][0], (int)points[i   ][1], 
			(int)points[nexI][0], (int)points[nexI][1] );
	}




}