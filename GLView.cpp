#include "StdAfx.h"
#include "GLView.h"


CGLView::CGLView(void)
{
	m_pGLDC = NULL;
}


CGLView::~CGLView(void)
{
}


// �������ƺ���
void CGLView::RenderScene(COpenGLDC* pDC)
{
	pDC->DrawCoord();
}


void CGLView::OnDraw(CDC* pDC)
{
	// TODO: �ڴ����ר�ô����/����û���
	if (m_pGLDC)
	{
		m_pGLDC->Ready();
		RenderScene(m_pGLDC);
		m_pGLDC->Finish();
	}
}


BOOL CGLView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: �ڴ����ר�ô����/����û���
	cs.style |= WS_CLIPSIBLINGS|WS_CLIPCHILDREN;

	return CView::PreCreateWindow(cs);
}

IMPLEMENT_DYNCREATE(CGLView, CView)

BEGIN_MESSAGE_MAP(CGLView, CView)
	ON_WM_CREATE()
	ON_WM_DESTROY()
	ON_WM_SIZE()
	ON_WM_ERASEBKGND()
END_MESSAGE_MAP()


int CGLView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CView::OnCreate(lpCreateStruct) == -1)
		return -1;

	// TODO:  �ڴ������ר�õĴ�������
	m_pGLDC = new COpenGLDC(this->GetSafeHwnd());//����COpenGLDC�Ķ���
	m_pGLDC->InitDC();                           //��ʼ������

	return 0;
}


void CGLView::OnDestroy()
{
	CView::OnDestroy();

	// TODO: �ڴ˴������Ϣ����������
	if(m_pGLDC) delete m_pGLDC;                   //�ͷŶ���
}


void CGLView::OnSize(UINT nType, int cx, int cy)
{
	CView::OnSize(nType, cx, cy);

	// TODO: �ڴ˴������Ϣ����������
	if(m_pGLDC) 
		m_pGLDC->GLResize(cx,cy);
}


BOOL CGLView::OnEraseBkgnd(CDC* pDC)
{
	// TODO: �ڴ������Ϣ�����������/�����Ĭ��ֵ

	//return CView::OnEraseBkgnd(pDC);
	return TRUE;
}


// ��ȡ��ǰģ�͵������ݺ�
BOOL CGLView::GetBox(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1)
{
	return FALSE;
}


// ���ų���,������1Ϊ��С��С��1Ϊ�Ŵ�
void CGLView::Zoom(double dScale)
{
	m_pGLDC->m_Camera.zoom(dScale);
	Invalidate();                 //ˢ����ͼ
}


// ����һ�����ʵ����űȣ��Խ�ģ��ȫ����ʾ�ڳ�����
void CGLView::ZoomAll()
{
	
}


// ʹ�õ����ӽ����۲�ģ��
void CGLView::OnViewType(UINT type)
{
	ASSERT(type >= VIEW_FRONT&&type <= VIEW_NW_ISOMETRIC);
	m_pGLDC->m_Camera.set_view_type(type);
	Invalidate();
}


// ����ǰ�����ߴ�İٷֱ��ƶ�����������dpx��dpy�ķ�Χ��0~1
void CGLView::MoveView(double dpx, double dpy)
{
	m_pGLDC->m_Camera.move_view(dpx,dpy);
	Invalidate();
}
