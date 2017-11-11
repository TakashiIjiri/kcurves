#include "stdafx.h"
#include "KCurveUI.h"
#include "MainForm.h"


using namespace kcurves;

KCurveUI::KCurveUI()
{
	m_bR = m_bL = m_bM = false;

}




void KCurveUI::LBtnDown (const int x, const int y)
{
	m_bL = true;
	
	m_dragPtId = pickCP(x,y);
	if(m_dragPtId == -1) m_CPs.push_back( EVec2d(x,y) );

	compute_kCurves( m_CPs, 10, m_kCurveCP, m_curves );
	MainForm_repaint();
}

void KCurveUI::LBtnUp   (const int x, const int y)
{
	m_dragPtId = -1;
	m_bL = false;
}




void KCurveUI::RBtnDown (const int x, const int y){
	m_bR = true;

}

void KCurveUI::MBtnDown (const int x, const int y){
	m_bM = true;
}

void KCurveUI::RBtnUp   (const int x, const int y){
	m_bR = false;
}


void KCurveUI::MBtnUp   (const int x, const int y){
	m_bM = false;
}

void KCurveUI::MouseMove(const int x, const int y){
	if( !m_bR && !m_bL && !m_bM ) return;

	if( m_dragPtId != -1) m_CPs[m_dragPtId] << x,y;
	compute_kCurves( m_CPs, 20, m_kCurveCP, m_curves );
	MainForm_repaint();	
}