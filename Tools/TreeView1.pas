unit TreeView1;

interface

uses
  System.SysUtils, System.Classes, Vcl.Controls, Vcl.ComCtrls,generaltypes;

type

  TTreeView1 = class(TTreeView)
  private
  Fcurve:Tcurve;
    { Déclarations privées }
  protected
    { Déclarations protégées }
  public
    { Déclarations publiques }
    procedure Setcurve(const Value:Tcurve);
  published
  property Curve:Tcurve read Fcurve write SetCurve;
    { Déclarations publiées }
  end;

procedure Register;

implementation

procedure Register;
begin
  RegisterComponents('Pairings', [TTreeView1]);
end;

{ TTreeView1 }

procedure TTreeView1.Setcurve(const Value: Tcurve);
begin
Fcurve:=Value;
value.GenerateParamsTreeView(TTreeView(Self));
end;

end.
