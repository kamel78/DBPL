unit Unit2;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, GeneralTypes, BNCurves,LargeIntegers, Vcl.StdCtrls;

type
  TForm2 = class(TForm)
    Button1: TButton;
    Edit1: TEdit;
    Edit2: TEdit;
    BNCurvePairing1: TBNCurvePairing;
    Edit3: TEdit;
    Edit4: TEdit;
    Label1: TLabel;
    Label2: TLabel;
    Edit5: TEdit;
    Label3: TLabel;
    Label4: TLabel;
    Edit6: TEdit;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Edit7: TEdit;
    Label8: TLabel;
    procedure Button1Click(Sender: TObject);
  private
    { Déclarations privées }
  public
    { Déclarations publiques }
  end;

var
  Form2: TForm2;

implementation

{$R *.dfm}

procedure TForm2.Button1Click(Sender: TObject);
var p,ap:G1BN;
    q,bq:G2BN;
    e,e1,e2:GTBN;
    a,b:lint;
    tim:ttime;
    h,m,s,ms:word;
begin
tim:=now;
p:=BNCurvePairing1.HashToG1Point('kamel');
q:=BNCurvePairing1.HashToG2Point('ali');
e:=BNCurvePairing1.Paire(p,q);
edit7.Text:=e.ToHexString;
GetRandomLIntLowerThan(a,BNCurvePairing1.CurveParams.r);
edit5.Text:=a.ToHexString;
GetRandomLIntLowerThan(b,BNCurvePairing1.CurveParams.r);
edit6.Text:=b.ToHexString;
ap.SetCurveParams(BNCurvePairing1.CurveParams);
bq.SetCurveParams(BNCurvePairing1.CurveParams);
ap:=a*p;
bq:=b*q;

e1:=BNCurvePairing1.Paire(ap,bq);
edit2.Text:=e1.ToHexString;
tim:=now-tim;
DecodeTime(tim,h,m,s,ms);
e2:=e.Pow(a*b);
edit1.Text:=e2.ToHexString;
label8.Caption:='Time :'+'(in '+inttostr(m)+':'+inttostr(s)+':'+inttostr(ms)+'ms)';
end;

end.
