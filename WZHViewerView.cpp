#include "stdafx.h"
#ifndef SHARED_HANDLERS
#include "WZHViewer.h"
#endif

#include "WZHViewerDoc.h"
#include "WZHViewerView.h"
#include "GLView.h"
#include "GeomBase.h"
#include "solid.h"

#include <propkey.h>
#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CstlView

IMPLEMENT_DYNCREATE(CWZHViewerView, CGLView)  /////////////

	BEGIN_MESSAGE_MAP(CWZHViewerView, CGLView)   ///////////////
		// 标准打印命令
		ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
		ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
		ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CWZHViewerView::OnFilePrintPreview)
		ON_WM_CONTEXTMENU()
		ON_WM_RBUTTONUP()
		ON_COMMAND(ID_stlin, &CWZHViewerView::Onstlin)
		ON_COMMAND(ID_ZOOMALL, &CWZHViewerView::OnZoomall)
		ON_COMMAND(ID_FRONT, &CWZHViewerView::OnFront)
		ON_COMMAND(ID_BACK, &CWZHViewerView::OnBack)
		ON_COMMAND(ID_LEFT, &CWZHViewerView::OnLeft)
		ON_COMMAND(ID_RIGHT, &CWZHViewerView::OnRight)
		ON_COMMAND(ID_TOP, &CWZHViewerView::OnTop)
		ON_COMMAND(ID_BOTTOM, &CWZHViewerView::OnBottom)
		ON_COMMAND(ID_SW, &CWZHViewerView::OnSw)
		ON_COMMAND(ID_SE, &CWZHViewerView::OnSe)
		ON_COMMAND(ID_NW, &CWZHViewerView::OnNw)
		ON_COMMAND(ID_NE, &CWZHViewerView::OnNe)
		ON_COMMAND(ID_SHADOW, &CWZHViewerView::OnShadow)
		ON_WM_KEYDOWN()
//		ON_WM_MOUSEHWHEEL()
ON_WM_MOUSEWHEEL()
ON_WM_LBUTTONDOWN()
ON_WM_LBUTTONUP()
ON_WM_MOUSEMOVE()
ON_WM_RBUTTONDOWN()
ON_COMMAND(ID_creatpoint, &CWZHViewerView::Oncreatpoint)
ON_COMMAND(ID_erazePoint, &CWZHViewerView::Onerazepoint)
ON_COMMAND(ID_getTXT, &CWZHViewerView::Ongettxt)
ON_COMMAND(ID_AdaBoostTXT, &CWZHViewerView::OnAdaboosttxt)
ON_COMMAND(ID_newfindPOINT, &CWZHViewerView::Onnewfindpoint)
ON_COMMAND(ID_Skeletonizing, &CWZHViewerView::OnSkeletonizing)
ON_COMMAND(ID_clear, &CWZHViewerView::Onclear)
ON_COMMAND(ID_SkeletonizingClose, &CWZHViewerView::OnSkeletonizingclose)
ON_COMMAND(ID_Skeletonizing_smoothing, &CWZHViewerView::OnSkeletonizingsmoothing)
ON_COMMAND(ID_sortSmooth, &CWZHViewerView::Onsortsmooth)
ON_COMMAND(ID_drawClose, &CWZHViewerView::Ondrawclose)
ON_COMMAND(ID_SDF, &CWZHViewerView::OnSdf)
ON_COMMAND(ID_SDF_GMM, &CWZHViewerView::OnSdfGmm)
ON_COMMAND(ID_DIANYUN, &CWZHViewerView::OnDianyun)
ON_COMMAND(ID_LBACK, &CWZHViewerView::OnLback)
ON_COMMAND(ID_directFIND, &CWZHViewerView::Ondirectfind)
ON_COMMAND(ID_initPath, &CWZHViewerView::Oninitpath)
	END_MESSAGE_MAP()

	// CstlView 构造/析构

	CWZHViewerView::CWZHViewerView()
	{
		// TODO: 在此处添加构造代码
		pSTLModel = 0;
		double x0=y0=z0=x1=y1=z1;
		m_LeftButtonDown = FALSE;
		m_OnRButtonDown=FALSE;
		m_xRotate = 0.0;
		m_yRotate = 0.0;
		m_point=0;
		m_erazePoint=0;
		Creat=0;
		Erase=0;
		Creat_DOWN=0;
		sortsmooth=0;
		Skeletonizingsmoothing=0;
		Lback=1;
		Lchange=0;
	}

	CWZHViewerView::~CWZHViewerView()
	{
		if (pSTLModel != 0)
		{
			delete pSTLModel;
			pSTLModel = 0;
		}
	}

	BOOL CWZHViewerView::PreCreateWindow(CREATESTRUCT& cs)
	{
		// TODO: 在此处通过修改
		//  CREATESTRUCT cs 来修改窗口类或样式

		return CView::PreCreateWindow(cs);
	}

	// CstlView 绘制

	//void CstlView::OnDraw(CDC* /*pDC*/)
	//{
	//	CstlDoc* pDoc = GetDocument();
	//	ASSERT_VALID(pDoc);
	//	if (!pDoc)
	//		return;
	//
	//	// TODO: 在此处为本机数据添加绘制代码
	//}


	/*CstlView 打印*/

	void CWZHViewerView::OnFilePrintPreview()
	{
#ifndef SHARED_HANDLERS
		AFXPrintPreview(this);
#endif
	}

	BOOL CWZHViewerView::OnPreparePrinting(CPrintInfo* pInfo)
	{
		// 默认准备
		return DoPreparePrinting(pInfo);
	}

	void CWZHViewerView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
	{
		// TODO: 添加额外的打印前进行的初始化过程
	}

	void CWZHViewerView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
	{
		// TODO: 添加打印后进行的清理过程
	}

	void CWZHViewerView::OnRButtonUp(UINT /* nFlags */, CPoint point)
	{
		ClientToScreen(&point);
		OnContextMenu(this, point);
	}

	void CWZHViewerView::OnContextMenu(CWnd* /* pWnd */, CPoint point)
	{
#ifndef SHARED_HANDLERS
		theApp.GetContextMenuManager()->ShowPopupMenu(IDR_POPUP_EDIT, point.x, point.y, this, TRUE);
#endif
	}


	// CstlView 诊断

#ifdef _DEBUG
	void CWZHViewerView::AssertValid() const
	{
		CView::AssertValid();
	}

	void CWZHViewerView::Dump(CDumpContext& dc) const
	{
		CView::Dump(dc);
	}

	CWZHViewerDoc* CWZHViewerView::GetDocument() const // 非调试版本是内联的
	{
		ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CWZHViewerDoc)));
		return (CWZHViewerDoc*)m_pDocument;
	}
#endif //_DEBUG


	// CstlView 消息处理程序


	
	void CWZHViewerView::RenderScene(COpenGLDC* pDC)
	{
		glRotatef(m_yRotate,0.0f, 1.0f, 0.0f); // Rock Z
		glRotatef(m_xRotate,1.0f, 0.0f, 0.0f); // Roll X

		glPushMatrix();
		CWZHViewerDoc* pDoc = GetDocument();
		ASSERT(pDoc);
		pDC->DrawCoord();
		//pDC->RenderInformation(pSTLModel->m_vecPFacetTri.size(),pSTLModel->m_vecPVert.size());
		if (pSTLModel != 0) ////////////////////////////////
		{
		
			if (!pSTLModel->IsEmpty()) //////////////////////////////
			{
				
					if (m_point==1)
					{
						pSTLModel->Draw2(pDC);
					}
					if(Skeletonizingsmoothing==1)		
					{
						//pSTLModel->Draw_obb(pDC,vec_obb);
						pSTLModel->Draw_ccPath(pDC,m_vecPPolyPolygons);
						pSTLModel->Draw1(pDC,Rbox);
					}
					 if(SDF==1)
					{
						pSTLModel->Draw3(pDC,Rbox);
					}
					if (Skeletonizingsmoothing==0)
					 pSTLModel->Draw(pDC);
			}
		}

	}
	
	

	void CWZHViewerView::Onstlin()
	{
		CFileDialog dlg(TRUE, "stl", NULL, OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT, 
		"Stereo Lithographic File(*.stl)|*.stl||",NULL);

	CString strFile;
	if (dlg.DoModal()==IDOK)
	  {
		pSTLModel = new CSTLModel();
		LPCTSTR file;
		file=dlg.GetPathName();

		strFile = dlg.GetPathName();
		pSTLModel->LoadSTLFile(file);
		OnZoomall();
		if(pSTLModel->IsEmpty())
			delete pSTLModel;
	   }
	}


   void CWZHViewerView::OnZoomall()
	{
	   pSTLModel->GetBox(x0,y0,z0,x1,y1,z1);
		m_pGLDC->m_Camera.zoom_all(pSTLModel->m_one.x0,pSTLModel->m_one.y0,pSTLModel->m_one.z0,pSTLModel->m_one.x1,pSTLModel->m_one.y1,pSTLModel->m_one.z1);
		Invalidate();
	}


   void CWZHViewerView::OnFront()
   {
	   initRotate();
	   OnViewType(VIEW_FRONT);// TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::OnBack()
   {
	   initRotate();
	   OnViewType(VIEW_BACK);// TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::OnLeft()
   {
	   initRotate();
	   OnViewType(VIEW_LEFT);// TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::OnRight()
   {
	   initRotate();
	  OnViewType(VIEW_RIGHT); // TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::OnTop()
   {
	   initRotate();
	   OnViewType(VIEW_TOP);// TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::OnBottom()
   {
	   initRotate();
	   OnViewType(VIEW_BOTTOM);// TODO: 在此添加命令处理程序代码

   }



   void CWZHViewerView::OnSw()
   {
	   initRotate();
	   OnViewType(VIEW_SW_ISOMETRIC);// TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::OnSe()
   {
	   initRotate();
	  OnViewType(VIEW_SE_ISOMETRIC); // TODO: 在此添加命令处理程序代码
   }



   void CWZHViewerView::OnNw()
   {
	   initRotate();
	   OnViewType(VIEW_NW_ISOMETRIC);// TODO: 在此添加命令处理程序代码
   }



   void CWZHViewerView::OnNe()
   {
	   initRotate();
	   OnViewType(VIEW_NE_ISOMETRIC);// TODO: 在此添加命令处理程序代码
   }



   void CWZHViewerView::OnShadow()
   {
	   m_pGLDC->Shading(!m_pGLDC->IsShading());
	   Invalidate(); // TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
   {
	   switch(nChar){
	   case VK_UP:
		   MoveView(0.0,0.02);
		   break;
	   case VK_DOWN:
		   MoveView(0.0,-0.02);
		   break;
	   case VK_RIGHT:
		   MoveView(0.02,0);
		   break;
	   case VK_LEFT:
		   MoveView(-0.02,0);
		   break;
	   }// TODO: 在此添加消息处理程序代码和/或调用默认值

	   CGLView::OnKeyDown(nChar, nRepCnt, nFlags);
   }


//   void CWZHViewerView::OnMouseHWheel(UINT nFlags, short zDelta, CPoint pt)
//   {
//	   // 此功能要求 Windows Vista 或更高版本。
//	   // _WIN32_WINNT 符号必须 >= 0x0600。
//	   // TODO: 在此添加消息处理程序代码和/或调用默认值
//	   if(zDelta > 100){
//		   Zoom(0.9);
//		   Invalidate();
//	   }
//	   if(zDelta < -100){
//		   Zoom(1.1);
//		   Invalidate();}
//
//	   CGLView::OnMouseHWheel(nFlags, zDelta, pt);
//   }


   BOOL CWZHViewerView::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)
   {
	   // TODO: 在此添加消息处理程序代码和/或调用默认值
	   if(zDelta > 100){
		   Zoom(0.9);
		   Invalidate();
	   }
	   if(zDelta < -100){
		   Zoom(1.1);
		   Invalidate();}

	   return CGLView::OnMouseWheel(nFlags, zDelta, pt);
   }


   void CWZHViewerView::OnLButtonDown(UINT nFlags, CPoint point)
   {
	   // TODO: 在此添加消息处理程序代码和/或调用默认值
	   m_LeftButtonDown = TRUE;
	   m_LeftDownPos = point;	

	   CGLView::OnLButtonDown(nFlags, point);
   }


   void CWZHViewerView::OnLButtonUp(UINT nFlags, CPoint point)
   {
	   // TODO: 在此添加消息处理程序代码和/或调用默认值
	   m_LeftButtonDown = FALSE;
	   m_erazePoint=0;
	  Creat=0;
	   Erase=0;
	   CGLView::OnLButtonUp(nFlags, point);
   }



   void CWZHViewerView::OnMouseMove(UINT nFlags, CPoint Apoint)
   {
	   // TODO: 在此添加消息处理程序代码和/或调用默认值
	   if(m_LeftButtonDown)
	   {
		   CSize rotate = m_LeftDownPos - Apoint;//计算出鼠标按下的时刻与现在，鼠标坐标点的差
	   m_LeftDownPos = Apoint;
	   //if(nFlags==MK_CONTROL)
	   if (m_erazePoint==0)
	   {
		   if(nFlags & MK_CONTROL)
			   MoveView(-0.001*rotate.cx,0.001*rotate.cy);

		   else
		   {m_yRotate -= rotate.cx;
		   m_xRotate -= rotate.cy;}
		   InvalidateRect(NULL,FALSE);
	   }
	   if (m_erazePoint==1)
	   {
		   POINT3D PP;
		   PP=ScreenToPoint(Apoint);
		   point vi=point(PP.x,PP.y,PP.z);
		   double juli=pCalcubase->dis(*(pSTLModel->m_vecPFacetTri.back()->m_PVerts[0]),*(pSTLModel->m_vecPFacetTri.back()->m_PVerts[1]));
		   const float *Match=KD->closest_to_pt(vi,juli);
		   int my_ID=0;
		   if (Match)
		   {
			   my_ID=(Match-&vec_p[0][0])/3;
		   }
		  if (Erase==1&&Creat==0)
		  {
			  pSTLModel->m_vecPVert[my_ID]->bused=0;
		  }
		  if (Erase==0&&Creat==1)
		  {
			  if (pSTLModel->m_vecPVert[my_ID]->bused==0)
			  {
				 /* if (pSTLModel->dis((*pSTLModel->m_vecPVert[my_ID]),(*templine.back()))<10*juli)
				  {*/
				     templine.push_back(pSTLModel->m_vecPVert[my_ID]);
				 // }
			  }
			 pSTLModel->m_vecPVert[my_ID]->bused=1;
		  }
		  Invalidate();
	   }
	   }	

	   CGLView::OnMouseMove(nFlags, Apoint);
   }
   void CWZHViewerView::initRotate()
   {
	   m_yRotate = 0;
	   m_xRotate = 0;
   }


   void CWZHViewerView::OnRButtonDown(UINT nFlags, CPoint point)
   {
	   // TODO: 在此添加消息处理程序代码和/或调用默认值
	   m_OnRButtonDown=TRUE;
	   ScreenToPoint(point);

   }
   POINT3D CWZHViewerView::ScreenToPoint(CPoint P)
   {
	   int hits;
	   UINT items[64];
	   double x,y,z;
	   x=y=z=0;
	   POINT3D PP;
	   PP=m_pGLDC->BeginSelection(P.x, P.y);
	   //RenderScene(m_pGLDC);
	   hits = m_pGLDC->EndSelection(items);
	   /*m_pGLDC->m_Camera.GetWxyz(x,y,z);
	   wp=CPoint3D(x,y,z);*/
	   /*wp.x=x;
	   wp.y=y;
	   wp.z=z;*/
	   return PP;
   }





   void CWZHViewerView::Oncreatpoint()
   {
	   vec_p.clear();
	   //pSTLModel->findEDGE();
	   Creat=1;
	   m_erazePoint=1;
	   Creat_DOWN=1;
	   m_point=1;
	   templine.clear();
	   for (int i=0;i<pSTLModel->m_vecPVert.size();i++)
	   {
		   point vi(pSTLModel->m_vecPVert[i]->x,pSTLModel->m_vecPVert[i]->y,pSTLModel->m_vecPVert[i]->z);
		   vec_p.push_back(vi);
	   }
	   float* v0=&vec_p[0][0];
	   KD = new KDtree(v0,pSTLModel->m_vecPVert.size());
	   Invalidate();
	   // TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::Onerazepoint()
   {
	   Erase=1;
	   m_erazePoint=1;
	   m_point=1;
	   Invalidate();
	   // TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::Ongettxt()
   {
	   CTXT pTXT;
	   pTXT.outputTXT(pSTLModel);
	   // TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::OnAdaboosttxt()
   {
	   CTXT pTXT;
	   pTXT.inputTXT(pSTLModel);
	   m_point = 1;
	   Invalidate();
	   // TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::Onnewfindpoint()
   {
	   // TODO: 在此添加命令处理程序代码
	    CWZHViewerDoc* pDoc = GetDocument();
		DoInput theINput;
		theINput.DoModal();
		theINput.theLIMIT;
	   if (Lchange==0)
	   {
		   pDoc->m_vecLp.clear();
		   pDoc->pCLaplacian->change(pSTLModel->m_vecPVert,pSTLModel,pDoc->m_vecLp,pDoc->m_vecV);
		   pDoc->pCboundary->featurPOINT(pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPHEdge,pSTLModel->m_vecPVert,theINput.theLIMIT);
		   pDoc->pCLaplacian->back(pSTLModel->m_vecPVert,pSTLModel,pDoc->m_vecLp,pDoc->m_vecV);
		   pDoc->pCboundary->featurPOINT_new(pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPHEdge,pSTLModel->m_vecPVert);
		   m_point=1;
		   Invalidate();
	   }
	  if (Lchange==1)
	  {
		  pDoc->pCboundary->featurPOINT(pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPHEdge,pSTLModel->m_vecPVert,theINput.theLIMIT);
		  m_point=1;
		  Invalidate();
	  }
	  
   }


   void CWZHViewerView::OnSkeletonizing()
   {
	   // TODO: 在此添加命令处理程序代码
	   CWZHViewerDoc* pDoc = GetDocument();
	   pDoc->pCboundary->Skeletonizing(RBLOOP,pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPHEdge,pSTLModel->m_vecPVert,pSTLModel->m_vecPEdge);
	   m_point=1;
	   Invalidate(); 
   }

   void CWZHViewerView::Onclear()
   {
	   // TODO: 在此添加命令处理程序代码
	   m_point=0;
	   m_erazePoint=0;
	   sortsmooth=0;
	   RBLOOP.clear();
	   for (int i=0;i<pSTLModel->m_vecPVert.size();i++)
	   {
		   pSTLModel->m_vecPVert[i]->bused=0;
		   pSTLModel->m_vecPVert[i]->bStatus=0;
	   }
	   for (int i=0;i<pSTLModel->m_vecPHEdge.size();i++)
	   {
		   pSTLModel->m_vecPHEdge[i]->bStatus=0;
		   pSTLModel->m_vecPHEdge[i]->bused=0;
	   }
	   templine.clear();
	   Invalidate();
   }


   void CWZHViewerView::OnSkeletonizingclose()
   {
	   // TODO: 在此添加命令处理程序代码
	   CWZHViewerDoc* pDoc = GetDocument();
	   DoInput2 theINPUT;
	   theINPUT.DoModal();
	   pDoc->pCboundary->Skeletonizing_close(RBLOOP,pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPVert,pSTLModel->m_vecPHEdge,theINPUT.theANG); 
		m_point=1;
		Invalidate();
   }


   void CWZHViewerView::OnSkeletonizingsmoothing()
   {
	   // TODO: 在此添加命令处理程序代码
	 CWZHViewerDoc* pDoc = GetDocument();
	 ///////////////////计时
	/* time_t start,ends;
	 start=clock();
	 pDoc->m_vecLp.clear();
	 pDoc->pCLaplacian->change(pSTLModel->m_vecPVert,pSTLModel,pDoc->m_vecLp,pDoc->m_vecV);
	 double INput=2.6;
	 pDoc->pCboundary->featurPOINT(pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPHEdge,pSTLModel->m_vecPVert,INput);
	 pDoc->pCLaplacian->back(pSTLModel->m_vecPVert,pSTLModel,pDoc->m_vecLp,pDoc->m_vecV);
	 pDoc->pCboundary->featurPOINT_new(pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPHEdge,pSTLModel->m_vecPVert);
	 OnSkeletonizing();
	 INput=60;
	 pDoc->pCboundary->Skeletonizing_close(RBLOOP,pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPVert,pSTLModel->m_vecPHEdge,INput); */
	 ///////////////////////////////////////////
	 vector<vector<bool>> IDlist;
	 //////////////////////////
	 pDoc->pCboundary->Skeletonizing_smoothing(Rbox,pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPVert,RBLOOP,pSTLModel->m_vecPHEdge,IDlist);
	 delete pSTLModel;
	 pSTLModel = new CSTLModel();
	 pSTLModel->LoadSTLFile("e:\\A.stl");
	 for (int i = 0; i < pSTLModel->m_vecPFacetTri.size();i++)
	 {
		 for (int k = 0; k < 3;++k)
		 {
			 pSTLModel->m_vecPFacetTri[i]->m_PVerts[k]->bused = IDlist[i][k];
		 }		 
	 }
	 pDoc->pCboundary->resetLOOP(RBLOOP, pSTLModel->m_vecPVert);
	 
	 ///////////////////////////////
	 pDoc->pCboundary->Skeletonizing_partlize(Rbox, pSTLModel->m_vecPFacetTri, pSTLModel->m_vecPVert);	 
	 for (int i=0;i<Rbox.size();i++)
	 {
		
		 OBB ob;
		 COBB::makeOBB(Rbox[i],ob);
		 ob.ID = i;
		 vec_obb.push_back(ob);
	 }
	 vec_PFACETTRI vecpFac;
	 for (int i = 0; i<RBLOOP.size(); i++)
	 {
		 Ccalcubase::FindOneRFac(RBLOOP[i][RBLOOP[i].size()%2],vecpFac);
		 for (int j=0;j<vecpFac.size();j++)
		 {
			int aaa=vecpFac[j]->ID;
			int bbb=vecpFac[(j+1)%vecpFac.size()]->ID;
			 if (aaa!=bbb)
			 {
				 if (vec_obb[aaa].isBranch==TRUE)
				 {
					 vec_obb[aaa].pBloop=&RBLOOP[i];
					 vec_obb[aaa].preID = bbb;
				 }
				 if (vec_obb[bbb].isBranch==TRUE)
				 {
					 vec_obb[bbb].pBloop=&RBLOOP[i];
					 vec_obb[bbb].preID = aaa;
				 }
				 break;
			 }
		 }
	 }
		   m_point=1;
		   Skeletonizingsmoothing=1;
		   Invalidate();
		   ///////////////////////////
		 /*ends=clock();
		   double time1=difftime(ends,start);
		   time1=0.001*time1;
		   CString s;
		   s.Format(_T("%f"),time1);
		   AfxMessageBox(s);*/
		   ////////////////////////////
   }


   void CWZHViewerView::Onsortsmooth()
   {
	   // TODO: 在此添加命令处理程序代码
	  // pSTLModel->sort_smooth(templine);
	    CWZHViewerDoc* pDoc = GetDocument();
		pDoc->pCboundary->sort_smooth(templine);
	   m_point=1;
	   sortsmooth=1;
	   //templine.clear();
	   Invalidate();
   }


   void CWZHViewerView::Ondrawclose()
   {
	    CWZHViewerDoc* pDoc = GetDocument();
	   if (sortsmooth==1)
	   {
		   Creat_DOWN=1;
		   DoInput2 theINPUT;
		   theINPUT.DoModal();
		   pDoc->pCboundary->Skeletonizing_link(templine,pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPVert,pSTLModel->m_vecPHEdge,theINPUT.theANG);
		  // pDoc->pCboundary->snaking(templine,pSTLModel->m_vecPVert);
		  /* for (int i=0;i<templine.size();i++)
		   {
			   if (templine[i]==templine[(i+1)%templine.size()])
			   {
				   templine.erase(templine.begin()+i);
				   i--;
			   }
		   }
		   pDoc->pCboundary->snake_smooth_close(templine);*/
		   RBLOOP.push_back(templine);
		   templine.clear();
		   Invalidate();
	   }
	   else
	   {
		   pDoc->pCboundary->sort_smooth(templine);
		   m_point=1;
		   DoInput2 theINPUT;
		   theINPUT.DoModal();
		   pDoc->pCboundary->Skeletonizing_link(templine,pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPVert,pSTLModel->m_vecPHEdge,theINPUT.theANG);
		   pDoc->pCboundary->snaking(templine,pSTLModel->m_vecPVert);
		//   pDoc->pCboundary->snake_smooth(templine);
		 //  pDoc->pCboundary->Skeletonizing_link_new(templine,pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPVert);
		   templine.clear();
		   Invalidate();
	   }
	   // TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::OnSdf()
   {
	   // TODO: 在此添加命令处理程序代码
	    CWZHViewerDoc* pDoc = GetDocument();
		pDoc->pCsdf->SDF_kmeans(8,RBLOOP,pSTLModel->m_vecPFacetTri,pSTLModel->vec_box,Rbox);
	   m_point=1;
	   Skeletonizingsmoothing=1;
	   Invalidate();
   }


   void CWZHViewerView::OnSdfGmm()
   {
	   CWZHViewerDoc* pDoc = GetDocument();
	  // pSTLModel->UpdateBox();
	   pDoc->pCsdf->SDF_GMM(3,RBLOOP,pSTLModel->m_vecPFacetTri,pSTLModel->vec_box,Rbox);
	   //pSTLModel->SDF_GMM(3,RBLOOP);
	   m_point=1;
	   Skeletonizingsmoothing=1;
	   Invalidate();
	   // TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::OnDianyun()
   {
	   CWZHViewerDoc* pDoc = GetDocument();
	  /* pDoc->pCsdf->SDF_getPOINT(pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPVert,pSTLModel->vec_box);
	   m_point=1;
	   Skeletonizingsmoothing=1;*/ 
	   Lchange=1;
	   Lback=0;
	   pDoc->m_vecLp.clear();
	   pDoc->pCLaplacian->change(pSTLModel->m_vecPVert,pSTLModel,pDoc->m_vecLp,pDoc->m_vecV);
	   Invalidate();
	   // TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::OnLback()
   {
	   Lchange=0;
	   Lback=1;
	    CWZHViewerDoc* pDoc = GetDocument();
		pDoc->pCLaplacian->back(pSTLModel->m_vecPVert,pSTLModel,pDoc->m_vecLp,pDoc->m_vecV);
		Invalidate();
	   // TODO: 在此添加命令处理程序代码
   }


   void CWZHViewerView::Ondirectfind()
   {
	   CWZHViewerDoc* pDoc = GetDocument();
	   DoInput theINput;
	   theINput.DoModal();
	   theINput.theLIMIT;
	   // TODO: 在此添加命令处理程序代码
	   pDoc->pCboundary->featurPOINT(pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPHEdge,pSTLModel->m_vecPVert,theINput.theLIMIT);
	   /////////////////////计时
	   /*time_t start,ends;
	   start=clock();
	   double INput=2.6;
	  pDoc->pCboundary->featurPOINT(pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPHEdge,pSTLModel->m_vecPVert,INput);
	  OnSkeletonizing();
	  INput=60;
	  pDoc->pCboundary->Skeletonizing_close(RBLOOP,pSTLModel->m_vecPFacetTri,pSTLModel->m_vecPVert,pSTLModel->m_vecPHEdge,INput); 
	  m_point=1;
	  Skeletonizingsmoothing=1;
	  ends=clock();
	  double time1=difftime(ends,start);
	  time1=0.001*time1;
	  CString s;
	  s.Format(_T("%f"),time1);
	  AfxMessageBox(s);*/
	   ///////////////////////
	   m_point=1;
	   Invalidate();
   }


   void CWZHViewerView::Oninitpath()
   {
	   // TODO: 在此添加命令处理程序代码
	   CWZHViewerDoc* pDoc = GetDocument();
	   for (vec_PFACETTRI vec_fac:Rbox)
	   {
		   for (PFACETTRI Fac:vec_fac)
		   {
			   Fac->becut = 0;
		   }
	   }
	   for (int i=0;i<vec_obb.size();i++)
	   {
		   if (vec_obb[i].isBranch==TRUE)
		   {
			   double r(0.05),h(0.1);
			   pDoc->pCcutterPath->get_CCpath(vec_obb[i],m_vecPPolyPolygons,r,h,Rbox);
		   } 
	   }
	   Invalidate();
   }
