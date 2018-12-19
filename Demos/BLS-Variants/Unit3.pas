unit Unit3;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls, Vcl.ComCtrls,
  System.ImageList, Vcl.ImgList;

type
  TForm3 = class(TForm)
    Button1: TButton;
    Button2: TButton;
    imgToolbar: TImageList;
    Panel2: TPanel;
    Panel3: TPanel;
    Label2: TLabel;
    ComboBox2: TComboBox;
    Button6: TButton;
    ComboBox1: TComboBox;
    Panel1: TPanel;
    GroupBox1: TGroupBox;
    RadioButton1: TRadioButton;
    RadioButton2: TRadioButton;
    RadioButton3: TRadioButton;
    RichEdit1: TRichEdit;
    procedure WndProc(var Msgs: TMessage); override;
    procedure Split(const Delimiter: Char; Input: string; const Strings: TStrings) ;
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure FormShow(Sender: TObject);
    procedure RadioButton1Click(Sender: TObject);
    procedure Button6Click(Sender: TObject);
  private
    { Déclarations privées }
  public
    { Déclarations publiques }
  end;

var
  Form3: TForm3;

implementation
 uses  shellapi, richedit, Unit1, Unit2;
{$R *.dfm}

procedure TForm3.Button1Click(Sender: TObject);
begin
ModalResult:=1;
end;

procedure TForm3.Button2Click(Sender: TObject);
begin
modalresult:=2;
end;

procedure TForm3.Button6Click(Sender: TObject);
begin
case ComboBox2.ItemIndex of
0:form1.TreeView11.Curve:=Form2.BNCurvePairing1;
1:form1.TreeView11.Curve:=Form2.KSS18CurvePairing1;
2:form1.TreeView11.Curve:=Form2.BLS24CurvePairing1;
end;
Form1.ShowModal;
end;

procedure TForm3.FormCreate(Sender: TObject);
var
  mask : Integer;
begin
  mask := SendMessage(RichEdit1.Handle, EM_GETEVENTMASK, 0, 0);
  SendMessage(RichEdit1.Handle, EM_SETEVENTMASK, 0, mask or ENM_LINK);
  SendMessage(RichEdit1.Handle, EM_AUTOURLDETECT, Integer(True), 0);
end;

procedure TForm3.FormShow(Sender: TObject);
begin
RadioButton1Click(self);
end;

procedure TForm3.RadioButton1Click(Sender: TObject);
var i:integer;
begin
if  RadioButton1.Checked then i:=0
else if RadioButton2.Checked then i:=1
else if RadioButton3.Checked then i:=2;
ComboBox1.ItemIndex:=i;
case i of
0:RichEdit1.Text:='Boneh, Dan, Ben Lynn, and Hovav Shacham. "Short signatures from the Weil pairing." In International Conference on the Theory and Application of Cryptology and Information Security, pp. 514-532. '+'Springer, Berlin, Heidelberg, 2001. https://link.springer.com/article/10.1007/s00145-004-0314-9';
1:RichEdit1.Text:='Boneh, Dan, and Xavier Boyen. "Short signatures without random oracles." In International Conference on the Theory and Applications of Cryptographic Techniques, pp. 56-73. Springer, Berlin,'+' Heidelberg, 2004. https://link.springer.com/chapter/10.1007/978-3-540-24676-3_4';
2:RichEdit1.Text:='Boneh, Dan, and Xavier Boyen. "Short signatures without random oracles." In International Conference on the Theory and Applications of Cryptographic Techniques, pp. 56-73. Springer, Berlin,'+' Heidelberg, 2004. https://link.springer.com/chapter/10.1007/978-3-540-24676-3_4';
end;
end;

procedure TForm3.Split(const Delimiter: Char; Input: string; const Strings: TStrings) ;
begin
   Assert(Assigned(Strings)) ;
   Strings.Clear;
   Strings.Delimiter := Delimiter;
   Strings.DelimitedText := Input;
end;
procedure Tform3.WndProc(var Msgs: TMessage);
var
  p: TENLink;
  strURL, link, path : string;
  A : TStringList;
begin
  if (Msgs.Msg = WM_NOTIFY) then
  begin
    if (PNMHDR(Msgs.LParam).code = EN_LINK) then
    begin
      p := TENLink(Pointer(TWMNotify(Msgs).NMHdr)^);
      if (p.msg = WM_LBUTTONDOWN) then
      begin
        SendMessage(RichEdit1.Handle, EM_EXSETSEL, 0, LongInt(@(p.chrg)));
        link := RichEdit1.SelText;
        A := TStringList.Create;
        Split(':', link, A) ;
        path := GetEnvironmentVariable('TEMP');
        strURL := 'file:\\' + path + '\' + A[1] + '.png';
        A.Free;
        ShowMessage(strURL);
        ShellExecute(Handle, 'open', PChar(strURL), 0, 0, SW_SHOWNORMAL);
      end
    end
  end;

  inherited;
end;


end.
