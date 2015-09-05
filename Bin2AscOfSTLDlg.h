// Bin2AscOfSTLDlg.h : header file
//

#if !defined(AFX_BIN2ASCOFSTLDLG_H__B7F20DD8_8034_4EBB_BE72_577495EAC76E__INCLUDED_)
#define AFX_BIN2ASCOFSTLDLG_H__B7F20DD8_8034_4EBB_BE72_577495EAC76E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

/////////////////////////////////////////////////////////////////////////////
// CBin2AscOfSTLDlg dialog

class CBin2AscOfSTLDlg : public CDialog
{
// Construction
public:
	CBin2AscOfSTLDlg(CWnd* pParent = NULL);	// standard constructor
	void OnButton1();
// Dialog Data
	//{{AFX_DATA(CBin2AscOfSTLDlg)
	enum { IDD = IDD_BIN2ASCOFSTL_DIALOG };
	CListBox	m_lst;
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CBin2AscOfSTLDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	//{{AFX_MSG(CBin2AscOfSTLDlg)
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	//afx_msg void OnButton1();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_BIN2ASCOFSTLDLG_H__B7F20DD8_8034_4EBB_BE72_577495EAC76E__INCLUDED_)
