#pragma once
#include "afxwin.h"

#include "OpenGLDC.h"

class CGLView :public CView
{
protected:
	COpenGLDC*  m_pGLDC;

protected:
	CGLView(void);
	DECLARE_DYNCREATE(CGLView)

public:
	virtual ~CGLView(void);

	// 场景绘制函数
	virtual void RenderScene(COpenGLDC* pDC);

public:
	virtual void OnDraw(CDC* pDC);    //overridden to draw this view
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);

protected:
	DECLARE_MESSAGE_MAP()
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnDestroy();
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);

	
	virtual BOOL GetBox(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1);// 获取当前模型的最大包容盒
	void Zoom(double dScale);     // 缩放场景
	void ZoomAll();           // 计算一个合适的缩放比，以将模型全部显示在场景中
	void OnViewType(UINT type);   // 使用典型视角来观察模型
	void MoveView(double dpx, double dpy);// 按当前场景尺寸的百分比移动场景，参数dpx、dpy的范围是0~1
};

