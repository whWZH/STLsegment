
// WZHViewerView.h : CWZHViewerView ��Ľӿ�
//

#include "gl\GL.h"
#include "gl\glaux.h"
#include "gl\GLU.h"
#include "gl\glut.h"
#include "GLView.h"
#include "KD-tree/KDtree.h"
#include "KD-tree/Vec.h"
#include "TXT.h";
#include "calcubase.h"
#include "DoInput.h"
#include "DoInput2.h"


#pragma once


class CWZHViewerView : public CGLView
{
protected: // �������л�����
	CWZHViewerView();
	DECLARE_DYNCREATE(CWZHViewerView)
	//CSTLModel* m_pSTLModel;

	// ����
public:
	CWZHViewerDoc* GetDocument() const;
	CSTLModel* pSTLModel;          ////////////////////////////////////////////
	Ccalcubase* pCalcubase;
	KDtree* KD;
	vector<point> vec_p;
	double x0,y0,z0,x1,y1,z1;
	BOOL m_LeftButtonDown;
	BOOL m_OnRButtonDown;
	CPoint m_LeftDownPos;
	float m_xRotate;
	float m_yRotate;
	BOOL m_point;//��ʾ������
	BOOL m_erazePoint;//ȥ������������
	BOOL Skeletonizingsmoothing;//�����ɫ
	BOOL Creat;
	BOOL Erase;
	BOOL Creat_DOWN;//����ȫ������
	BOOL sortsmooth;
	BOOL SDF;//SDF����

	BOOL Lchange;
	BOOL Lback;

	POINT3D ScreenToPoint(CPoint P);
	vector<vector<PVERT>> RBLOOP;
	vector<vector<PFACETTRI>> RBfac;
	vector<vector<PFACETTRI>> Rbox;
	vec_PVERT templine;
	vector<OBB> vec_obb;

	vec_PPOLYPOLYGON	m_vecPPolyPolygons;
	//double getz;//�ָ�ƽ��
	// ����
public:

	// ��д
public:
	//virtual void OnDraw(CDC* pDC);  // ��д�Ի��Ƹ���ͼ   ////////////////////////////////
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);

	// ʵ��
public:
	virtual ~CWZHViewerView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

	// ���ɵ���Ϣӳ�亯��
protected:
	afx_msg void OnFilePrintPreview();
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	DECLARE_MESSAGE_MAP()
public:
	virtual void RenderScene(COpenGLDC* pDC);
	afx_msg void Onstlin();
	afx_msg void OnZoomall();
	afx_msg void OnFront();
	afx_msg void OnBack();
	afx_msg void OnLeft();
	afx_msg void OnRight();
	afx_msg void OnTop();
	afx_msg void OnBottom();
	afx_msg void OnSw();
	afx_msg void OnSe();
	afx_msg void OnNw();
	afx_msg void OnNe();
	afx_msg void OnShadow();
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
//	afx_msg void OnMouseHWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
private:
	void initRotate();
public:
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void Oncreatpoint();
	afx_msg void Onerazepoint();
	afx_msg void Ongettxt();
	afx_msg void OnAdaboosttxt();
	afx_msg void Onnewfindpoint();
	afx_msg void OnSkeletonizing();
	afx_msg void Onclear();
	afx_msg void OnSkeletonizingclose();
	afx_msg void OnSkeletonizingsmoothing();
	afx_msg void Onsortsmooth();
	afx_msg void Ondrawclose();
	afx_msg void OnSdf();
	afx_msg void OnSdfGmm();
	afx_msg void OnDianyun();
	afx_msg void OnLback();
	afx_msg void Ondirectfind();
	afx_msg void Oninitpath();
};

#ifndef _DEBUG  // stlView.cpp �еĵ��԰汾
inline CWZHViewerDoc* CWZHViewerView::GetDocument() const
{ return reinterpret_cast<CWZHViewerDoc*>(m_pDocument); }
#endif

