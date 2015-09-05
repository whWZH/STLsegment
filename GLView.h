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

	// �������ƺ���
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

	
	virtual BOOL GetBox(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1);// ��ȡ��ǰģ�͵������ݺ�
	void Zoom(double dScale);     // ���ų���
	void ZoomAll();           // ����һ�����ʵ����űȣ��Խ�ģ��ȫ����ʾ�ڳ�����
	void OnViewType(UINT type);   // ʹ�õ����ӽ����۲�ģ��
	void MoveView(double dpx, double dpy);// ����ǰ�����ߴ�İٷֱ��ƶ�����������dpx��dpy�ķ�Χ��0~1
};

