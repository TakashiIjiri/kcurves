#include "stdafx.h"
#include "KCurveUI.h"
#include "MainForm.h"
#include "catmullromspline.h"


using namespace kcurves;

KCurveUI::KCurveUI()
{
	m_bR = m_bL = m_bM = false;
}


void KCurveUI::UpdateCurve()
{
	if (MainForm_isClosed())
		compute_kCurves(m_cps, 15, m_kcurve_cps, m_curves);
	else
		compute_kCurves_open(m_cps, 15, m_kcurve_cps, m_curves);

	m_catmullrom_curve = catmullromspline::compute_catmullrom_spline(m_cps, 0.5f,30);
}



void KCurveUI::LBtnDown (const int x, const int y)
{
	m_bL = true;
	
	m_dragPtId = pickCP(x,y);
	if(m_dragPtId == -1) m_cps.push_back( EVec2f(x,y) );
	UpdateCurve();
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

	if( m_dragPtId != -1) m_cps[m_dragPtId] << (float)x, (float)y;
	UpdateCurve();
	MainForm_repaint();	
}