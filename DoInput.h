#pragma once


// DoInput 对话框

class DoInput : public CDialogEx
{
	DECLARE_DYNAMIC(DoInput)

public:
	DoInput(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~DoInput();

// 对话框数据
	enum { IDD = IDD_DIALOG1 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()
public:
	double theLIMIT;
};
