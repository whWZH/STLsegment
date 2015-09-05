
// WZHViewerDoc.h : CWZHViewerDoc 类的接口
//


#include "STLModel.h"
#include "Boundary.h"
#include "SDF.h"
#include "Laplacian.h"
#include "CutterPath.h"

#pragma once


class CWZHViewerDoc : public CDocument
{
protected: // 仅从序列化创建
	CWZHViewerDoc();
	DECLARE_DYNCREATE(CWZHViewerDoc)

	// 特性
public:
	CSTLModel* pSTLModel;//////////////////////////////////////////
	CSDF* pCsdf;
	Cboundary* pCboundary;
	CLaplacian* pCLaplacian;
	CcutterPath* pCcutterPath;
	vector<POINT3D> m_vecLp;
	vec_VECTOR3D m_vecV;
	// 操作
public:
	

	// 重
public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
#ifdef SHARED_HANDLERS
	virtual void InitializeSearchContent();
	virtual void OnDrawThumbnail(CDC& dc, LPRECT lprcBounds);
#endif // SHARED_HANDLERS

	// 实现
public:
	virtual ~CWZHViewerDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

	// 生成的消息映射函数
protected:
	DECLARE_MESSAGE_MAP()

#ifdef SHARED_HANDLERS
	// 用于为搜索处理程序设置搜索内容的 Helper 函数
	void SetSearchContent(const CString& value);
#endif // SHARED_HANDLERS
public:
	afx_msg void OnBack();
};
