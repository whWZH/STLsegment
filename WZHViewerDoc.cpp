
// WZHViewerDoc.cpp : CWZHViewerDoc ���ʵ��
//

#include "stdafx.h"
// SHARED_HANDLERS ������ʵ��Ԥ��������ͼ������ɸѡ�������
// ATL ��Ŀ�н��ж��壬�����������Ŀ�����ĵ����롣
#ifndef SHARED_HANDLERS
#include "WZHViewer.h"
#endif

#include "WZHViewerDoc.h"

#include <propkey.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CstlDoc

IMPLEMENT_DYNCREATE(CWZHViewerDoc, CDocument)

	BEGIN_MESSAGE_MAP(CWZHViewerDoc, CDocument)
	END_MESSAGE_MAP()


	// CstlDoc ����/����

	CWZHViewerDoc::CWZHViewerDoc()
	{
		pSTLModel = 0;
		// TODO: �ڴ����һ���Թ������

	}

	CWZHViewerDoc::~CWZHViewerDoc()
	{
		if (pSTLModel != 0)
		{
			delete pSTLModel;
			pSTLModel = 0;
		}
	}

	BOOL CWZHViewerDoc::OnNewDocument()
	{
		if (!CDocument::OnNewDocument())
			return FALSE;

		// TODO: �ڴ�������³�ʼ������
		// (SDI �ĵ������ø��ĵ�)

		return TRUE;
	}




	// CstlDoc ���л�

	void CWZHViewerDoc::Serialize(CArchive& ar)
	{
		if (ar.IsStoring())
		{
			// TODO: �ڴ���Ӵ洢����
		}
		else
		{
			// TODO: �ڴ���Ӽ��ش���
		}
	}

#ifdef SHARED_HANDLERS

	// ����ͼ��֧��
	void CstlDoc::OnDrawThumbnail(CDC& dc, LPRECT lprcBounds)
	{
		// �޸Ĵ˴����Ի����ĵ�����
		dc.FillSolidRect(lprcBounds, RGB(255, 255, 255));

		CString strText = _T("TODO: implement thumbnail drawing here");
		LOGFONT lf;

		CFont* pDefaultGUIFont = CFont::FromHandle((HFONT) GetStockObject(DEFAULT_GUI_FONT));
		pDefaultGUIFont->GetLogFont(&lf);
		lf.lfHeight = 36;

		CFont fontDraw;
		fontDraw.CreateFontIndirect(&lf);

		CFont* pOldFont = dc.SelectObject(&fontDraw);
		dc.DrawText(strText, lprcBounds, DT_CENTER | DT_WORDBREAK);
		dc.SelectObject(pOldFont);
	}

	// ������������֧��
	void CstlDoc::InitializeSearchContent()
	{
		CString strSearchContent;
		// ���ĵ����������������ݡ�
		// ���ݲ���Ӧ�ɡ�;���ָ�

		// ����:  strSearchContent = _T("point;rectangle;circle;ole object;")��
		SetSearchContent(strSearchContent);
	}

	void CstlDoc::SetSearchContent(const CString& value)
	{
		if (value.IsEmpty())
		{
			RemoveChunk(PKEY_Search_Contents.fmtid, PKEY_Search_Contents.pid);
		}
		else
		{
			CMFCFilterChunkValueImpl *pChunk = NULL;
			ATLTRY(pChunk = new CMFCFilterChunkValueImpl);
			if (pChunk != NULL)
			{
				pChunk->SetTextValue(PKEY_Search_Contents, value, CHUNK_TEXT);
				SetChunkValue(pChunk);
			}
		}
	}

#endif // SHARED_HANDLERS

	// CstlDoc ���

#ifdef _DEBUG
	void CWZHViewerDoc::AssertValid() const
	{
		CDocument::AssertValid();
	}

	void CWZHViewerDoc::Dump(CDumpContext& dc) const
	{
		CDocument::Dump(dc);
	}
#endif //_DEBUG


	// CstlDoc ����


