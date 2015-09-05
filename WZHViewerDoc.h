
// WZHViewerDoc.h : CWZHViewerDoc ��Ľӿ�
//


#include "STLModel.h"
#include "Boundary.h"
#include "SDF.h"
#include "Laplacian.h"
#include "CutterPath.h"

#pragma once


class CWZHViewerDoc : public CDocument
{
protected: // �������л�����
	CWZHViewerDoc();
	DECLARE_DYNCREATE(CWZHViewerDoc)

	// ����
public:
	CSTLModel* pSTLModel;//////////////////////////////////////////
	CSDF* pCsdf;
	Cboundary* pCboundary;
	CLaplacian* pCLaplacian;
	CcutterPath* pCcutterPath;
	vector<POINT3D> m_vecLp;
	vec_VECTOR3D m_vecV;
	// ����
public:
	

	// ��
public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
#ifdef SHARED_HANDLERS
	virtual void InitializeSearchContent();
	virtual void OnDrawThumbnail(CDC& dc, LPRECT lprcBounds);
#endif // SHARED_HANDLERS

	// ʵ��
public:
	virtual ~CWZHViewerDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

	// ���ɵ���Ϣӳ�亯��
protected:
	DECLARE_MESSAGE_MAP()

#ifdef SHARED_HANDLERS
	// ����Ϊ����������������������ݵ� Helper ����
	void SetSearchContent(const CString& value);
#endif // SHARED_HANDLERS
public:
	afx_msg void OnBack();
};
