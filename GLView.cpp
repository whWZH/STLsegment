#include "StdAfx.h"
#include "GLView.h"


CGLView::CGLView(void)
{
	m_pGLDC = NULL;
}


CGLView::~CGLView(void)
{
}


// 场景绘制函数
void CGLView::RenderScene(COpenGLDC* pDC)
{
	pDC->DrawCoord();
}


void CGLView::OnDraw(CDC* pDC)
{
	// TODO: 在此添加专用代码和/或调用基类
	if (m_pGLDC)
	{
		m_pGLDC->Ready();
		RenderScene(m_pGLDC);
		m_pGLDC->Finish();
	}
}


BOOL CGLView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: 在此添加专用代码和/或调用基类
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

	// TODO:  在此添加您专用的创建代码
	m_pGLDC = new COpenGLDC(this->GetSafeHwnd());//创建COpenGLDC的对象
	m_pGLDC->InitDC();                           //初始化对象

	return 0;
}


void CGLView::OnDestroy()
{
	CView::OnDestroy();

	// TODO: 在此处添加消息处理程序代码
	if(m_pGLDC) delete m_pGLDC;                   //释放对象
}


void CGLView::OnSize(UINT nType, int cx, int cy)
{
	CView::OnSize(nType, cx, cy);

	// TODO: 在此处添加消息处理程序代码
	if(m_pGLDC) 
		m_pGLDC->GLResize(cx,cy);
}


BOOL CGLView::OnEraseBkgnd(CDC* pDC)
{
	// TODO: 在此添加消息处理程序代码和/或调用默认值

	//return CView::OnEraseBkgnd(pDC);
	return TRUE;
}


// 获取当前模型的最大包容盒
BOOL CGLView::GetBox(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1)
{
	return FALSE;
}


// 缩放场景,（大于1为缩小，小于1为放大）
void CGLView::Zoom(double dScale)
{
	m_pGLDC->m_Camera.zoom(dScale);
	Invalidate();                 //刷新视图
}


// 计算一个合适的缩放比，以将模型全部显示在场景中
void CGLView::ZoomAll()
{
	
}


// 使用典型视角来观察模型
void CGLView::OnViewType(UINT type)
{
	ASSERT(type >= VIEW_FRONT&&type <= VIEW_NW_ISOMETRIC);
	m_pGLDC->m_Camera.set_view_type(type);
	Invalidate();
}


// 按当前场景尺寸的百分比移动场景，参数dpx、dpy的范围是0~1
void CGLView::MoveView(double dpx, double dpy)
{
	m_pGLDC->m_Camera.move_view(dpx,dpy);
	Invalidate();
}
