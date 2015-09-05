#pragma once


// DoInput2 对话框

class DoInput2 : public CDialogEx
{
	DECLARE_DYNAMIC(DoInput2)

public:
	DoInput2(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~DoInput2();

// 对话框数据
	enum { IDD = IDD_DIALOG2 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()
public:
	double theANG;
};
