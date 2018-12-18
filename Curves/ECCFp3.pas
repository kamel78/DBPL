unit ECCFp3;

interface

uses Fp3Arithmetic,System.SysUtils,LargeIntegers,Fp18Arithmetic,Fp6Arithmetic,HashFunctions, GeneralTypes,VCL.Dialogs;
 type

  Fp3Point=record     //  Definition of a Point on EC Over Fp3 (E(Fp3))
      public
      CurveParams: PtrCurveParams;
      X,Y,Z:Fp3Int;
      ComputeLigneValue:boolean;//KSS18
      ComputeLigneAtFp6:boolean;//MNT
      CurrentXp,CurrentYp:LInt;
      LineAtP:Fp18Int;
      LineAtQPFp6:Fp6Int;
      Infinity:boolean;
      procedure SetCurveParams(Value: PtrCurveParams;ComputLmd:boolean=false);
      procedure SetPointValue(Value:Fp3Int);
      function ToDecimalString:String;
      function ToHexString:String;
      procedure InitFromStrings(Xa,Xb,Xc,Ya,Yb,Yc:String;Za:String='1';Zb:String='0';Zc:String='0');
      procedure SetAsRandomPoint;
      procedure SetAsRandomTorsionPoint;
      procedure SetToDefaultGenerator;
      procedure SetAsTorsionFromHash(hash:TBytes);
      procedure SetAsTorsionFromString(s:String);
      procedure SetPairingPointCoordinates(PtX,PtY:LInt);
      function CompressToArray:TBytes;
      function IsOnTheCurve:boolean;
      Procedure DeCompressFromArray(a:TBytes);
      function FrobeniusMap_i(i:integer):Fp3Point;
      class operator Add(Left, Right: Fp3Point): Fp3Point;
      class operator Subtract(Left, Right: Fp3Point): Fp3Point;
      class operator Multiply(Left: LInt; Right: Fp3Point): Fp3Point;
      class operator Equal(const Left,Right: Fp3Point): Boolean;
      class operator NotEqual(const Left, Right: Fp3Point): Boolean;
      class operator Negative(const Value: Fp3Point): Fp3Point;
      end;

    procedure _Add_Jacobian_Fp3_Point(Left,Right:Fp3Point;var Result:Fp3Point;ComputeLignValue:boolean=false);
    procedure _Double_Jacobian_Fp3_Point(const Value:Fp3Point;var Result:Fp3Point;ComputeLignValue:boolean=false);
    procedure _Add_Projective_Fp3_Point(Left{Q},Right{T}:Fp3Point;var Result:Fp3Point;ComputeLignValue:boolean=false);
    procedure _Double_Projective_Fp3_Point(const Value:Fp3Point;var Result:Fp3Point;ComputeLigneValue:boolean=false);
    procedure _Add_Affine_Fp3_Point(Left,Right:Fp3Point;var Result:Fp3Point;ComputeLignValue:boolean=false);
    procedure _Double_Affine_Fp3_Point(const Value:Fp3Point;var Result:Fp3Point;ComputeLignValue:boolean=false);
    procedure _Double_Fp3_Point(Value:Fp3Point;var Result:Fp3Point;CoordSys:CoordinatesSystem=csAffine);
    procedure _Add_Fp3_Point(Left,Right:Fp3Point;var Result:Fp3Point;CoordSys:CoordinatesSystem=csAffine);
    procedure _Sub_Fp3_Point(Left,Right:Fp3Point;var Result:Fp3Point;CoordSys:CoordinatesSystem=csAffine);
    procedure _Evaluate_Fp3_Point(value:Fp3Int; var Result:Fp3Point);
    procedure _Neg_Fp3_Point(Value:Fp3Point;var Result:Fp3Point);
    function _Is_OnCurve_Fp3Point(value:Fp3Point):boolean;
    procedure _Frobenius_Map_i(const Value:Fp3Point;pow:integer;var Result:Fp3Point);
    procedure _Mul_Fp_Fp3Point(const Left:LInt;const Right:Fp3Point;var Result:Fp3Point;CoordSys:CoordinatesSystem=csAffine);
    procedure _Mul_NAF_Fp3_FpPoint(const Left:LInt;const Right:Fp3Point;var Result:Fp3Point;CoordSys:CoordinatesSystem=csAffine);
    procedure _Jacobian_To_Affine_Fp3Point(Value:Fp3Point;var Result:Fp3Point);
    procedure _Projective_To_Affine_Fp3Point(Value:Fp3Point;var Result:Fp3Point);
    function _Are_Equals_Fp3Points(left,right:Fp3Point):boolean;
    procedure _Fast_KSS18_Mul_Htw_Fp3Point(const Qx0:Fp3Point;var Result:Fp3Point);    //Right Should be different than Right

 implementation

{*******************************************************************************}
            ///      Procedures for Elliptic Curves Arithmetic on Fp3
{*******************************************************************************}

        {**********   Add two Fp3 Points  using Affine coordinates *****************}
procedure _Add_Affine_Fp3_Point(Left,Right:Fp3Point;var Result:Fp3Point;ComputeLignValue:boolean=false);
var t: array[0..4] of Fp3Int;
    tmp:Fp3Int;
begin
if Left.Infinity then Result:=Right
else if Right.Infinity then Result:=Left
else begin
     if _Compare_Fp3(Left.X,Right.X)=0 then begin
                            if _Compare_Fp3(Left.Y,Right.Y)=0 then _Double_Affine_Fp3_Point(Left,Result)
                            else begin
                                 Result.Z.a:=1;
                                 Result.Z.b:=0;
                                 Result.Z.c:=0;
                                 Result.SetCurveParams(Right.CurveParams,ComputeLignValue);
                                 Result.Infinity:=true;
                                 end;
                            end
     else begin
          Result.SetCurveParams(Left.CurveParams,ComputeLignValue);
          _Sub_Fp3(Right.Y,Left.Y,t[0]);
          _Sub_Fp3(Right.X,Left.X,t[1]);
          _Inv_Fp3(t[1],t[2]);
          _Mul_Fp3(t[0],t[2],t[4]); // Lambda
          _Sqr_Fp3(t[4],t[3]);
          if Result.ComputeLigneValue then begin    /// KSS18 Curves
          if Left.CurveParams.TwistMode=twDType then begin
                                                     Result.LineAtP.a.a.a:=Right.CurveParams.P -Right.CurrentYp;
                                                     Result.LineAtP.a.a.b:=0;
                                                     Result.LineAtP.a.a.c:=0;
                                                     Result.LineAtP.a.b.SetToZero;
                                                     Result.LineAtP.a.c.SetToZero;
                                                     _Mul_FP_Fp3(Right.CurrentXp,t[4],Result.LineAtP.b.a);
                                                     _Mul_Fp3(Left.x,t[4],t[0]);
                                                     _Sub_Fp3(Left.y,t[0],Result.LineAtP.b.b);
                                                     Result.LineAtP.b.c.SetToZero;
                                                     end
          else begin
                // KSS18 Curves
                     Result.LineAtP.b.b.a:=Right.CurveParams.P -Right.CurrentYp;
                     Result.LineAtP.b.b.b:=0;
                     Result.LineAtP.b.b.c:=0;
                     Result.LineAtP.b.a.SetToZero;
                     Result.LineAtP.a.c.SetToZero;
                     _Mul_FP_Fp3(Right.CurrentXp,t[4],Result.LineAtP.a.b);
                     _Mul_Fp3(Left.x,t[4],t[0]);
                     _Sub_Fp3(Left.y,t[0],Result.LineAtP.a.a);
                     Result.LineAtP.b.c.SetToZero;
               end;
          end
               else if Result.ComputeLigneAtFp6 then begin   //// MNT Curves
                                                 tmp.SetFieldParams(Result.CurveParams.FieldParam);
                                                 _Add_FP_Fp3(Right.CurrentXp,Left.X,tmp);
                                                 _Mul_Fp3(t[4],tmp,tmp);
                                                 _Sub_Fp3(Left.Y,tmp,tmp);
                                                 Result.LineAtQPFp6.a.a:=tmp.a;
                                                 Result.LineAtQPFp6.a.b:=Right.CurrentYp;

                                                 Result.LineAtQPFp6.b.a:=0;
                                                 Result.LineAtQPFp6.b.b:=tmp.c;
                                                 Result.LineAtQPFp6.c.a:=tmp.b;
                                                 Result.LineAtQPFp6.c.b:=0;

                                                 end;

          _Sub_Fp3(t[3],Left.X,Result.X);
          _Sub_Fp3(Result.X,Right.X,Result.X);
          _Sub_Fp3(Left.X,Result.X,t[2]);
          _Mul_Fp3(t[4],t[2],Result.Y);
          _Sub_Fp3(Result.Y,Left.Y,Result.Y);
          Result.Z.a:=1;
          Result.Z.b:=0;
          Result.Z.c:=0;
          end;
     end;
end;

        {**********   Double an Fp3 Point using Affine coordinates *****************}
procedure _Double_Affine_Fp3_Point(const Value:Fp3Point;var Result:Fp3Point;ComputeLignValue:boolean=false);
var t: array[0..3] of Fp3Int;
    ValX,ValY,tmp:Fp3Int;
begin
if(Value.Infinity)or(Value.Y.a.Data.i16[-2]=0)and (Value.Y.b.Data.i16[-2]=0) and (Value.Y.c.Data.i16[-2]=0) then Result.Infinity:=true
else begin
     Result.Infinity:=false;
     Result.SetCurveParams(Value.CurveParams,ComputeLignValue);
     ValX:=Value.X;
     ValY:=Value.Y;
     _Sqr_Fp3(Value.X,t[0]);
     _Mul_FP_Fp3(LInt(3),t[0],t[1]);
     _Add_Fp3(t[1],Value.CurveParams.AtwFp3,t[1]);
     _Mul_FP_Fp3(LInt(2),Value.Y,t[3]);
     _Inv_Fp3(t[3],t[3]);
     _Mul_Fp3(t[1],t[3],t[2]);  // Lambda
     _Sqr_Fp3(t[2],t[3]);
     _Sub_Fp3(t[3],Value.X,t[3]);
     _Sub_Fp3(t[3],Value.X,Result.X);
     _Sub_Fp3(ValX,Result.x,t[3]);
     _Mul_Fp3(t[2],t[3],Result.Y);
     _Sub_Fp3(Result.y,ValY,Result.Y);
     if Result.ComputeLigneValue then begin
      if Value.CurveParams.TwistMode=twDType then begin
                                                  _Sub_LInt(Value.CurveParams.P,Value.CurrentYp,Result.LineAtP.a.a.a);
                                                  Result.LineAtP.a.a.b:=0;
                                                  Result.LineAtP.a.a.c:=0;
                                                  Result.LineAtP.a.b.SetToZero;
                                                  Result.LineAtP.a.c.SetToZero;
                                                  _Mul_FP_Fp3(Value.CurrentXp,t[2],Result.LineAtP.b.a);
                                                  _Mul_Fp3(Valx,t[2],t[1]);
                                                  _Sub_Fp3(Valy,t[1],t[1]);
                                                  Result.LineAtP.b.b:=t[1];
                                                  Result.LineAtP.b.c.SetToZero;
                                                  end
      else begin
               _Sub_LInt(Value.CurveParams.P,Value.CurrentYp,Result.LineAtP.b.b.a);
               Result.LineAtP.b.b.b:=0;
               Result.LineAtP.b.b.c:=0;
               Result.LineAtP.b.a.SetToZero;
               Result.LineAtP.a.c.SetToZero;
               _Mul_FP_Fp3(Value.CurrentXp,t[2],Result.LineAtP.a.b);
               _Mul_Fp3(Valx,t[2],t[1]);
               _Sub_Fp3(Valy,t[1],t[1]);
               Result.LineAtP.a.a:=t[1];
               Result.LineAtP.b.c.SetToZero;

           end;
        end
      else           if Result.ComputeLigneAtFp6 then begin   //// MNT Curves
                                                 tmp.SetFieldParams(Result.CurveParams.FieldParam);
                                                 _Add_FP_Fp3(Value.CurrentXp,ValX,tmp);
                                                 _Mul_Fp3(t[2],tmp,tmp);
                                                 _Sub_Fp3(ValY,tmp,tmp);
                                                 Result.LineAtQPFp6.a.a:=tmp.a;
                                                 Result.LineAtQPFp6.a.b:=Value.CurrentYp;

                                                 Result.LineAtQPFp6.b.a:=0;
                                                 Result.LineAtQPFp6.b.b:=tmp.c;
                                                 Result.LineAtQPFp6.c.a:=tmp.b;
                                                 Result.LineAtQPFp6.c.b:=0;
                                            end

     end
end;
        {**********   Add two Fp3 Points using Jacobian coordinates *****************}
        // (X/Z^2,Y/Z^3)        https://eprint.iacr.org/2010/354.pdf
procedure _Add_Jacobian_Fp3_Point(Left,Right:Fp3Point;var Result:Fp3Point;ComputeLignValue:boolean=false);
 var t: array[0..10] of Fp3Int;
     Zr2,Ly2:Fp3Int;
begin
if Left.Infinity then Result:=Right
else if Right.Infinity then Result:=Left
else begin
     if _Compare_Fp3(Left.X,Right.X)=0 then begin
                            if _Compare_Fp3(Left.Y,Right.Y)=0 then
                            if _Compare_Fp3(Left.Z,Right.Z)=0 then _Double_Jacobian_Fp3_Point(Left,Result)
                            else begin
                                 Result.Z.a:=1;
                                 Result.Z.b:=0;
                                 Result.Z.c:=0;
                                 Result.SetCurveParams(Right.CurveParams,ComputeLignValue);
                                 Result.Infinity:=true;
                                 end;
                            end
     else begin
          Result.SetCurveParams(Left.CurveParams,ComputeLignValue);
          _Sqr_Fp3(Right.Z,Zr2);
          _Mul_Fp3(Left.X,Zr2,t[0]);
          _Add_Fp3(Left.Y,Right.Z,t[2]);
          _Sqr_Fp3(t[2],t[1]);
          _Sqr_Fp3(Left.Y,Ly2);
          _Sub_Fp3(t[1],Ly2,t[1]);
          _Sub_Fp3(t[1],Zr2,t[2]);
          _Mul_Fp3(t[2],Zr2,t[1]);
          _Sub_Fp3(t[0],Right.X,t[2]);
          _Sqr_Fp3(t[2],t[3]);
          _Mul_FP_Fp3(LInt(4),t[3],t[4]);
          _Mul_Fp3(t[4],t[2],t[5]);
          _Mul_FP_Fp3(LInt(2),Right.Y,t[6]);
          _Sub_Fp3(t[1],t[6],t[6]);
          _Mul_Fp3(t[6],Left.X,t[9]);
          _Mul_Fp3(Right.X,t[4],t[7]);
          _Mul_FP_Fp3(LInt(2),t[7],t[10]);
          _Sqr_Fp3(t[6],Result.X);
          _Sub_Fp3(Result.X,t[5],Result.x);
          _Sub_Fp3(Result.X,t[10],Result.X);
          _Add_Fp3(Right.Z,t[2],Result.Z);
          _Sqr_Fp3(Result.Z,t[10]);
          _Sub_Fp3(t[10],Zr2,Result.Z);
          _Sub_Fp3(Result.Z,t[3],Result.Z);
          _Add_Fp3(Left.Y,Result.Z,t[10]);
          _Sub_Fp3(t[7],Result.X,t[7]);
          _Mul_Fp3(t[7],t[6],t[8]);
          _Mul_Fp3(Right.Y,t[5],t[0]);
          _Shl_LInt(t[0].a,1);
          _Shl_LInt(t[0].b,1);
          _Shl_LInt(t[0].c,1);
          _Mod_LInt(t[0].a,Left.CurveParams.P,t[0].a);
          _Mod_LInt(t[0].b,Left.CurveParams.P,t[0].b);
          _Mod_LInt(t[0].c,Left.CurveParams.P,t[0].c);
          _Sub_Fp3(t[8],t[0],Result.Y);
          if Result.ComputeLigneValue then begin
          if Left.CurveParams.TwistMode=twDType then begin
                                                     _Sqr_Fp3(t[10],t[3]);
                                                     _Sub_Fp3(t[3],Ly2,t[3]);
                                                     _Sqr_Fp3(Result.Z,t[4]);
                                                     _Sub_Fp3(t[3],t[4],t[10]);
                                                     _Shl_LInt(t[9].a,1);
                                                     _Shl_LInt(t[9].b,1);
                                                     _Shl_LInt(t[9].c,1);
                                                     _Mod_LInt(t[9].a,Right.CurveParams.P,t[9].a);
                                                     _Mod_LInt(t[9].b,Right.CurveParams.P,t[9].b);
                                                     _Mod_LInt(t[9].c,Right.CurveParams.P,t[9].c);
                                                     _Sub_Fp3(t[9],t[10],t[9]);
                                                     _Mul_FP_Fp3(Right.CurrentYp,Result.Z,t[10]);
                                                     _Shl_LInt(t[10].a,1);
                                                     _Shl_LInt(t[10].b,1);
                                                     _Shl_LInt(t[10].c,1);
                                                     _Mod_LInt(t[10].a,Right.CurveParams.P,t[10].a);
                                                     _Mod_LInt(t[10].b,Right.CurveParams.P,t[10].b);
                                                     _Mod_LInt(t[10].c,Right.CurveParams.P,t[10].c);
                                                     _Sub_LInt(Right.CurveParams.P,t[6].a,t[6].a);
                                                     _Sub_LInt(Right.CurveParams.P,t[6].b,t[6].b);
                                                     _Sub_LInt(Right.CurveParams.P,t[6].c,t[6].c);
                                                     _Mul_FP_Fp3(Right.CurrentXp,t[6],t[1]);
                                                     _Shl_LInt(t[1].a,1);
                                                     _Shl_LInt(t[1].b,1);
                                                     _Shl_LInt(t[1].c,1);
                                                     _Mod_LInt(t[1].a,Right.CurveParams.P,t[1].a);
                                                     _Mod_LInt(t[1].b,Right.CurveParams.P,t[1].b);
                                                     _Mod_LInt(t[1].c,Right.CurveParams.P,t[1].c);
                                                     Result.LineAtP.a.a:=t[10];
                                                     Result.LineAtP.a.b.SetToZero;
                                                     Result.LineAtP.a.c.SetToZero;
                                                     Result.LineAtP.b.a:=t[1];
                                                     Result.LineAtP.b.b:=t[9];
                                                     Result.LineAtP.b.c.SetToZero;
                                                     end
          else begin
               _Sqr_Fp3(t[10],t[3]);
               _Sub_Fp3(t[3],Ly2,t[3]);
               _Sqr_Fp3(Result.Z,t[4]);
               _Sub_Fp3(t[3],t[4],t[10]);
               _Shl_LInt(t[9].a,1);
               _Shl_LInt(t[9].b,1);
               _Shl_LInt(t[9].c,1);
               _Mod_LInt(t[9].a,Right.CurveParams.P,t[9].a);
               _Mod_LInt(t[9].b,Right.CurveParams.P,t[9].b);
               _Mod_LInt(t[9].c,Right.CurveParams.P,t[9].c);
               _Sub_Fp3(t[9],t[10],t[9]);
               _Mul_FP_Fp3(Right.CurrentYp,Result.Z,t[10]);
               _Shl_LInt(t[10].a,1);
               _Shl_LInt(t[10].b,1);
               _Shl_LInt(t[10].c,1);
               _Mod_LInt(t[10].a,Right.CurveParams.P,t[10].a);
               _Mod_LInt(t[10].b,Right.CurveParams.P,t[10].b);
               _Mod_LInt(t[10].c,Right.CurveParams.P,t[10].c);
               _Sub_LInt(Right.CurveParams.P,t[6].a,t[6].a);
               _Sub_LInt(Right.CurveParams.P,t[6].b,t[6].b);
               _Sub_LInt(Right.CurveParams.P,t[6].c,t[6].c);
               _Mul_FP_Fp3(Right.CurrentXp,t[6],t[1]);
               _Shl_LInt(t[1].a,1);
               _Shl_LInt(t[1].b,1);
               _Shl_LInt(t[1].c,1);
               _Mod_LInt(t[1].a,Right.CurveParams.P,t[1].a);
               _Mod_LInt(t[1].b,Right.CurveParams.P,t[1].b);
               _Mod_LInt(t[1].c,Right.CurveParams.P,t[1].c);
               Result.LineAtP.b.b:=t[10];
               Result.LineAtP.b.a.SetToZero;
               Result.LineAtP.a.c.SetToZero;
               Result.LineAtP.a.b:=t[1];
               Result.LineAtP.a.a:=t[9];
               Result.LineAtP.b.c.SetToZero;
               end;
            end;
          end;
     end;
end;
        {**********   Double an Fp3 Point using Jacobian coordinates *****************}
        // (X/Z^2,Y/Z^3)        https://eprint.iacr.org/2010/354.pdf
procedure _Double_Jacobian_Fp3_Point(const Value:Fp3Point;var Result:Fp3Point;ComputeLignValue:boolean=false);
var t: array[0..7] of Fp3Int;
    ZQ2:Fp3Int;
begin
if(Value.Infinity)or((Value.Y.a.Data.i16[-2]=0) and (Value.Y.b.Data.i16[-2]=0) and (Value.Y.c.Data.i16[-2]=0)) then Result.Infinity:=true
else begin
     Result.Infinity:=false;
     Result.SetCurveParams(Value.CurveParams,ComputeLignValue);
     _Sqr_Fp3(Value.Z,ZQ2);
     _Sqr_Fp3(Value.X,t[0]);
     _Sqr_Fp3(Value.Y,t[1]);
     _Sqr_Fp3(t[1],t[2]);
     _Add_Fp3(t[1],Value.X,t[3]);
     _Sqr_Fp3(t[3],t[4]);
     _Sub_Fp3(t[4],t[0],t[3]);
     _Sub_Fp3(t[3],t[2],t[3]);
     _Shl_LInt(t[3].a,1);
     _Shl_LInt(t[3].b,1);
     _Shl_LInt(t[3].c,1);
     _Mod_LInt(t[3].a,Value.CurveParams.P,t[3].a);
     _Mod_LInt(t[3].b,Value.CurveParams.P,t[3].b);
     _Mod_LInt(t[3].c,Value.CurveParams.P,t[3].c);
     _Mul_FP_Fp3(LInt(3),t[0],t[4]);
     _Add_Fp3(Value.X,t[4],t[6]);
     _Sqr_Fp3(t[4],t[5]);
     _Mul_FP_Fp3(LInt(2),t[3],Result.X);
     _Sub_Fp3(t[5],Result.X,Result.X);
     _Add_Fp3(Value.Y,Value.Z,t[7]);
     _Sqr_Fp3(t[7],Result.Z);
     _Sub_Fp3(Result.Z,t[1],Result.Z);
     _Sub_Fp3(Result.Z,ZQ2,Result.Z);
     _Sub_Fp3(t[3],Result.X,Result.Y);
     _Mul_Fp3(Result.Y,t[4],t[7]);
     _Mul_FP_Fp3(LInt(8),t[2],Result.Y);
     _Sub_Fp3(t[7],Result.Y,Result.Y);
     if Result.ComputeLigneValue then begin
     if Value.CurveParams.TwistMode=twDType then begin
                                                 _Mul_Fp3(t[4],ZQ2,t[3]);
                                                 _Mul_FP_Fp3(LInt(2),t[3],t[7]);
                                                 _Sub_LInt(Value.CurveParams.P,t[7].a,t[7].a);
                                                 _Sub_LInt(Value.CurveParams.P,t[7].b,t[7].b);
                                                 _Sub_LInt(Value.CurveParams.P,t[7].c,t[7].c);
                                                 _Mul_FP_Fp3(Value.CurrentXp,t[7],Result.LineAtP.b.a);
                                                 _Sqr_Fp3(t[6],t[7]);
                                                 _Mul_FP_Fp3(LInt(4),t[1],t[2]);
                                                 _Sub_Fp3(t[7],t[0],t[7]);
                                                 _Sub_Fp3(t[7],t[5],t[7]);
                                                 _Sub_Fp3(t[7],t[2],Result.LineAtP.b.b);
                                                 _Mul_Fp3(Result.Z,ZQ2,t[2]);
                                                 _Mul_FP_Fp3(LInt(2),t[2],t[7]);
                                                 _Mul_FP_Fp3(Value.CurrentYp,t[7],Result.LineAtP.a.a);
                                                 Result.LineAtP.a.b.SetToZero;
                                                 Result.LineAtP.a.c.SetToZero;
                                                 Result.LineAtP.b.c.SetToZero;
                                                 end
     else begin
          _Mul_Fp3(t[4],ZQ2,t[3]);
          _Mul_FP_Fp3(LInt(2),t[3],t[7]);
          _Sub_LInt(Value.CurveParams.P,t[7].a,t[7].a);
          _Sub_LInt(Value.CurveParams.P,t[7].b,t[7].b);
          _Sub_LInt(Value.CurveParams.P,t[7].c,t[7].c);
          _Mul_FP_Fp3(Value.CurrentXp,t[7],Result.LineAtP.a.b);
          _Sqr_Fp3(t[6],t[7]);
          _Mul_FP_Fp3(LInt(4),t[1],t[2]);
          _Sub_Fp3(t[7],t[0],t[7]);
          _Sub_Fp3(t[7],t[5],t[7]);
          _Sub_Fp3(t[7],t[2],Result.LineAtP.a.a);
          _Mul_Fp3(Result.Z,ZQ2,t[2]);
          _Mul_FP_Fp3(LInt(2),t[2],t[7]);
          _Mul_FP_Fp3(Value.CurrentYp,t[7],Result.LineAtP.b.b);
          Result.LineAtP.b.a.SetToZero;
          Result.LineAtP.a.c.SetToZero;
          Result.LineAtP.b.c.SetToZero;
          end;
        end;
     end
end;

        {**********   Double an Fp3 Point (Projective coordinates (X/Z,Y/Z))*****************}
        // Costello et al's   https://cryptojedi.org/papers/edate-20100614.pdf
procedure _Double_Projective_Fp3_Point(const Value:Fp3Point;var Result:Fp3Point;ComputeLigneValue:boolean=false);
var tmp,tmp1,lam,A,B,C,D,E,F,G:Fp3Int;
begin
if(Value.Infinity)or((Value.Y.a.Data.i16[-2]=0)and (Value.Y.b.Data.i16[-2]=0))or(_Is_Fp3_Null(Value.Z)) then Result.Infinity:=true
else begin
     Result.Infinity:=false;
     Result.SetCurveParams(Value.CurveParams,ComputeLigneValue);
     _Sqr_Fp3(Value.X,A);
     _Sqr_Fp3(Value.Y,B);
     _Sqr_Fp3(Value.Z,C);
     _Mul_Fp3(C,Value.CurveParams.BtwFp3,tmp);
     _Nomod_Add_Fp3(tmp,tmp,D);
     _Add_Fp3(tmp,D,D);
     _Nomod_Add_Fp3(Value.X,Value.Y,tmp);
     _Sqr_Fp3(tmp,E);
     _Nomod_Sub_Fp3(E,A,E);
     _Sub_Fp3(E,B,E);
     _Nomod_Add_Fp3(Value.Y,Value.Z,tmp);
     _Sqr_Fp3(tmp,F);
     _Nomod_Sub_Fp3(F,B,F);
     _Sub_Fp3(F,C,F);
     _Nomod_Add_Fp3(D,D,G);
     _Add_Fp3(D,G,G);

     if Result.ComputeLigneValue then begin
     if Value.CurveParams.TwistMode=twDType then begin
                                                 _Mul_FP_Fp3(Value.CurrentYp,F,tmp);
                                                 _Neg_Fp3(tmp,Result.LineAtP.a.a);     //l11
                                                 Result.LineAtP.a.b.SetToZero;
                                                 Result.LineAtP.a.c.SetToZero;
                                                 _Nomod_Add_Fp3(A,A,tmp);
                                                 _Add_Fp3(A,tmp,tmp);
                                                 _Mul_FP_Fp3(Value.CurrentXp,tmp,Result.LineAtP.b.a);    //l02
                                                 _Sub_Fp3(D,B, Result.LineAtP.b.b);     //l00
                                                 Result.LineAtP.b.c.SetToZero;
                                                 end
     else begin
          _Mul_FP_Fp3(Value.CurrentYp,F,tmp);
          _Neg_Fp3(tmp,Result.LineAtP.b.b);     //l11
          Result.LineAtP.b.a.SetToZero;
          Result.LineAtP.a.c.SetToZero;
          _Nomod_Add_Fp3(A,A,tmp);
          _Add_Fp3(A,tmp,tmp);
          _Mul_FP_Fp3(Value.CurrentXp,tmp,Result.LineAtP.a.b);    //l02
          _Sub_Fp3(D,B, Result.LineAtP.a.a);     //l00
          Result.LineAtP.b.c.SetToZero;
         end;
      end
      else if Result.ComputeLigneAtFp6 then begin   /// MNT Curves
                                                /// Projective ccordinates for MNT Curves has not been implemented (For the Ate pairing)
                                                ///  since the same pointes are used for KSS18, and has not the same projective formulation
                                                ///  (the MNT curves has non nule A parameters) So implemented Projective curves for MNT in the same
                                                ///  unit will make the code very  harder to understand, In addition the gain with projective coordinates for MNT is not important (1 ms)
                                            end;
      _Nomod_Sub_Fp3(B,G,tmp);
     _Mul_Fp3(tmp,E,Result.X);
     _Add_Fp3(B,G,tmp);
     _Sqr_Fp3(tmp,Result.Y);
     _Sqr_Fp3(D,tmp);
     _Shl_LInt(tmp.a,2);
     _Shl_LInt(tmp.b,2);
     _Shl_LInt(tmp.c,2);
     _Sub_Fp3(Result.Y,tmp,Result.Y);
     _Nomod_Add_Fp3(tmp,tmp,tmp);
     _Sub_Fp3(Result.Y,tmp,Result.Y);
     _Mul_Fp3(B,F,Result.Z);
     _Shl_LInt(Result.Z.a,2);
     _Shl_LInt(Result.Z.b,2);
     _Shl_LInt(Result.Z.c,2);
     _Mod_LInt(Result.Z.a,Value.CurveParams.P,Result.Z.a);
     _Mod_LInt(Result.Z.b,Value.CurveParams.P,Result.Z.b);
     _Mod_LInt(Result.Z.c,Value.CurveParams.P,Result.Z.c);
   end
end;

        {**********   Add two Fp3 Points (Projective coordinates)*****************}
        // Mixed addition (Left.Z is always equal to 1)
        // Costello et al's   https://cryptojedi.org/papers/edate-20100614.pdf
procedure _Add_Projective_Fp3_Point(Left{Q},Right{T}:Fp3Point;var Result:Fp3Point;ComputeLignValue:boolean=false);
var A,B,C,D,tmp1,tmp2,lam:Fp3Int;
begin
if Left.Infinity then Result:=Right
else if Right.Infinity then Result:=Left
else begin
     if _Compare_Fp3(Left.X,Right.X)=0 then begin
                            if _Compare_Fp3(Left.Y,Right.Y)=0
                            then if _Compare_Fp3(Left.Z,Right.Z)=0 then _Double_Projective_Fp3_Point(Left,Result)
                            else begin
                                 Result.Infinity:=true;
                                 Result.Z.a:=1;
                                 Result.Z.b:=0;
                                 Result.Z.c:=0;
                                 Result.SetCurveParams(Right.CurveParams,ComputeLignValue);
                                 end;
                            end
     else begin
          Result.SetCurveParams(Left.CurveParams,ComputeLignValue);
          _Mul_Fp3(Right.Z,Left.X,A);
          _Sub_Fp3(Right.X,A,A);
          _Mul_Fp3(Right.Z,Left.Y,B);
          _Sub_Fp3(Right.Y,B,B);
          _Sqr_Fp3(A,tmp1);
          _Mul_Fp3(Right.X,tmp1,tmp2); // tmp2 used later as Right.X
          _Mul_Fp3(tmp1,A,C);
          _Mul_Fp3(B,B,tmp1);
          _Mul_Fp3(tmp1,Right.Z,D);
          _Nomod_Add_Fp3(D,C,D);
          _Add_Fp3(tmp2,tmp2,tmp1);
          _Sub_Fp3(D,tmp1,D);
          _Sub_Fp3(tmp2,D,tmp2);
          _Mul_Fp3(B,tmp2,Result.Y);
          _Mul_Fp3(Right.Y,C,tmp2);
          _Sub_Fp3(Result.Y,tmp2,Result.Y);
          _Mul_Fp3(A,D,Result.X);
          _Mul_Fp3(Right.Z,C,Result.Z);
          if Result.ComputeLigneValue then begin
          if  Left.CurveParams.TwistMode=twDType then begin
                                                           _Mul_FP_Fp3(-Right.CurrentYp,A,Result.LineAtP.a.a);      //L0,1
                                                           Result.LineAtP.a.b.SetToZero;
                                                           Result.LineAtP.a.c.SetToZero;
                                                           _Mul_FP_Fp3(Right.CurrentXp,B,Result.LineAtP.b.a);
                                                           _Mul_Fp3(B,Left.X,tmp1);
                                                           _Mul_Fp3(A,Left.Y,tmp2);
                                                           _Sub_Fp3(tmp2,tmp1,Result.LineAtP.b.b);
                                                           Result.LineAtP.b.c.SetToZero;
                                                           end
          else begin
               _Mul_FP_Fp3(-Right.CurrentYp,A,Result.LineAtP.b.b);      //L0,1
               Result.LineAtP.b.a.SetToZero;
               Result.LineAtP.a.c.SetToZero;
               _Mul_FP_Fp3(Right.CurrentXp,B,tmp1);
               _Mul_Fp3(B,Left.X,tmp1);
               _Mul_Fp3(A,Left.Y,tmp2);
               _Sub_Fp3(tmp2,tmp1,Result.LineAtP.a.a);
               Result.LineAtP.b.c.SetToZero;
               end;
            end
          else if Result.ComputeLigneAtFp6 then begin  /// MNT Curves
                                                /// Projective ccordinates for MNT Curves has not been implemented (For the Ate pairing)
                                                ///  since the same pointes are used for KSS18, and has not the same projective formulation
                                                ///  (the MNT curves has non nule A parameters) So implemented Projective curves for MNT in the same
                                                ///  unit will make the code very  harder to understand, In addition the gain with projective coordinates for MNT is not important (1 ms)
                                                end;
          end;
     end;
end;

        {**********   Add two Fp3 Points *****************}
procedure _Add_Fp3_Point(Left,Right:Fp3Point;var Result:Fp3Point;CoordSys:CoordinatesSystem=csAffine);
begin
case CoordSys of
  csAffine:_Add_Affine_Fp3_Point(Left,Right,Result);
  csProjective:_Add_Projective_Fp3_Point(Left,Right,Result);
  csJacobian:_Add_Jacobian_Fp3_Point(Left,Right,Result);
end;
end;

        {**********   Double two Fp3 Points *****************}
procedure _Double_Fp3_Point(Value:Fp3Point;var Result:Fp3Point;CoordSys:CoordinatesSystem=csAffine);
begin
case CoordSys of
  csAffine:_Double_Affine_Fp3_Point(Value,Result);
  csProjective:_Double_Projective_Fp3_Point(Value,Result);
  csJacobian:_Double_Jacobian_Fp3_Point(Value,Result);
end;
end;

        {**********   Substract two Fp3 Points *****************}
procedure _Sub_Fp3_Point(Left,Right:Fp3Point;var Result:Fp3Point;CoordSys:CoordinatesSystem=csAffine);
var tmp:Fp3Point;
begin
_neg_Fp3_point(Right,tmp);
case CoordSys of
  csAffine:_Add_Affine_Fp3_Point(Left,tmp,Result);
  csProjective:_Add_Projective_Fp3_Point(Left,tmp,Result);
  csJacobian:_Add_Jacobian_Fp3_Point(Left,tmp,Result);
end;
end;

       {**********   Convert Jacobian/Affine for an Fp3 Point *****}
        //(X/Z^2,Y/Z^3)
procedure _Jacobian_To_Affine_Fp3Point(Value:Fp3Point;var Result:Fp3Point);
var t1,t2:Fp3Int;
begin
Result.SetCurveParams(Value.CurveParams);
if (_Is_Fp3_Null(Value.Z)) then Result.Infinity:=true
else begin
     t1:=Value.z.Inverse;
     t2:=t1.Sqr;
     Result.X:=(Value.X*t2);
     Result.y:=(Value.y*(t1*t2));
     Result.z.a:=1;
     Result.z.b:=0;
     Result.z.c:=0;
     end;
end;

       {**********   Convert Projective/Affine for an Fp3 Point *****}
       // (X/Z,Y/Z)
procedure _Projective_To_Affine_Fp3Point(Value:Fp3Point;var Result:Fp3Point);
var t1,t2:Fp3Int;
begin
Result.SetCurveParams(Value.CurveParams);
if (_Is_Fp3_Null(Value.Z)) then Result.Infinity:=true
else begin
     t1:=Value.z.Inverse;
     Result.X:=(Value.X*t1);
     Result.y:=(Value.y*t1);
     Result.z.a:=1;
     Result.z.b:=0;
     Result.z.c:=0;
     end;
end;

       {**********   Normalize coordinates to Affine for an Fp3 Point *****}
procedure _Normalize_PF2_Point(value:Fp3Point;var Result:Fp3Point;CoordSys:CoordinatesSystem=csAffine);
begin
case  CoordSys of
  csAffine:result:=value;
  csJacobian:_Jacobian_To_Affine_Fp3Point(value,Result);
  csProjective:_Projective_To_Affine_Fp3Point(Value,Result);
end;
end;

      { **********   Multiply the Cofactor with an Fp3 Point ******* }
      // Fast Multiplication of A point by the Twist Cofactor (Generated By  Luis Dominguez)
      // For KSS Curves
      // https://eprint.iacr.org/2008/530.pdf
procedure _Fast_KSS18_Mul_Htw_Fp3Point(const Qx0:Fp3Point;var Result:Fp3Point);    //Right Should be different than Result
var t1,tmp,t2,t3,t4,t5,t6,Qx0_,Qx1_,Qx1,Qx2_,Qx2,Qx3:Fp3Point;
begin
Result.CurveParams :=Qx0.CurveParams;
if Qx0.CurveParams<>nil then  Result.CurveParams:=Qx0.CurveParams;
if Qx0.Infinity then begin
                     Result.Infinity:=true;
                     exit;
                     end;
Result.Infinity:=False;

_Neg_Fp3_Point(Qx0,Qx0_);
_Mul_Fp_Fp3Point(Qx0.CurveParams.u,Qx0,Qx1,csProjective);
_Neg_Fp3_Point(Qx1,Qx1_);
_Mul_Fp_Fp3Point(Qx0.CurveParams.u,Qx1,Qx2,csProjective);
_Neg_Fp3_Point(Qx2,Qx2_);
_Mul_Fp_Fp3Point(Qx0.CurveParams.u,Qx2,Qx3,csProjective);

t1:=Qx0;
_Frobenius_Map_i(Qx1_,2,t2);

_Frobenius_Map_i(Qx1,5,t3);
_Add_Affine_Fp3_Point(t3,Qx1,t3);

_Frobenius_Map_i(Qx2_,2,t4);
_Frobenius_Map_i(Qx2,1,tmp);
_Add_Affine_Fp3_Point(t4,tmp,t4);
_Frobenius_Map_i(Qx1,3,tmp);
_Add_Affine_Fp3_Point(t4,tmp,t4);
_Frobenius_Map_i(Qx0_,4,t5);

_Frobenius_Map_i(Qx3,1,t6);
_Frobenius_Map_i(Qx2,5,tmp);
_Add_Affine_Fp3_Point(t6,tmp,t6);
_Frobenius_Map_i(Qx2_,4,tmp);
_Add_Affine_Fp3_Point(t6,tmp,t6);
_Frobenius_Map_i(Qx0,3,tmp);
_Add_Affine_Fp3_Point(t6,tmp,t6);
_Frobenius_Map_i(Qx0,1,tmp);
_Add_Affine_Fp3_Point(t6,tmp,t6);

_Add_Affine_Fp3_Point(t2,t1,t2);
_Add_Affine_Fp3_Point(t1,t1,t1);
_Add_Affine_Fp3_Point(t1,t3,t1);
_Add_Affine_Fp3_Point(t1,t2,t1);
_Add_Affine_Fp3_Point(t4,t2,t4);
_Add_Affine_Fp3_Point(t5,t1,t5);
_Add_Affine_Fp3_Point(t4,t1,t4);
_Add_Affine_Fp3_Point(t5,t4,t5);
_Add_Affine_Fp3_Point(t4,t6,t4);
_Add_Affine_Fp3_Point(t5,t5,t5);
_Add_Affine_Fp3_Point(t4,t5,Result);
Result.Z.SetToOne;
end;

      { **********   Multiply the Cofactor with an Fp3 Point ******* }
      // Fast Multiplication of A point by the Twist Cofactor (Generated By  Luis Dominguez)
      // For MNT Curves
      // https://eprint.iacr.org/2008/530.pdf
procedure _Fast_MNT_Mul_Htw_Fp3Point(const Q:Fp3Point;var Result:Fp3Point);    //Right Should be different than Result
var tmp,tmp1,_2uQ:Fp3Point;
    u,us:Fp3Int;
    zz:Fp6Int;
begin
u.SetFieldParams(Q.CurveParams.FieldParam);
u.SetToZero;
u.b:=1;
Result.CurveParams:=Q.CurveParams;
if Q.CurveParams<>nil then  Result.CurveParams:=Q.CurveParams;
if Q.Infinity then begin
                     Result.Infinity:=true;
                     exit;
                     end;
Result.Infinity:=False;
tmp.SetCurveParams(Q.CurveParams);
_Double_Affine_Fp3_Point(Q,tmp);
_Mul_Fp_Fp3Point(Q.CurveParams.u,tmp,_2uQ,csAffine);
_Frobenius_Map_i(_2uQ,1,tmp);
_Frobenius_Map_i(_2uQ,2,tmp1);
_Add_Affine_Fp3_Point(tmp,tmp1,Result);
Result.Z.SetToOne;
end;


       {**********   Multiply a scalar with an Fp3 Point *******}
procedure _Mul_Fp_Fp3Point(const Left:LInt;const Right:Fp3Point;var Result:Fp3Point; CoordSys:CoordinatesSystem=csAffine);    //Right Should be different than Right
var i:integer;
begin
Result.CurveParams:=Right.CurveParams;
if (Right.Infinity)or (Left.Data.i16[-2]=0) then Result.Infinity:=true
else begin
     Result.Infinity:=true;
     if Left=2 then begin
                    case CoordSys of
                    csAffine:_Double_Affine_Fp3_Point(Right,Result);
                    csProjective:_Double_Projective_Fp3_Point(Right,Result);
                    csJacobian:_Double_Jacobian_Fp3_Point(Right,Result);
                    end;
                    end
     else
     for i:=Left.BitLength-1 downto 0 do begin
                                            case CoordSys of
                                              csAffine:_Double_Affine_Fp3_Point(Result,Result);
                                              csProjective:_Double_Projective_Fp3_Point(Result,Result);
                                              csJacobian:_Double_Jacobian_Fp3_Point(Result,Result);
                                            end;
                                         if _Is_BitSet_At(Left,i) then begin
                                                                  case CoordSys of
                                                                    csAffine:_Add_Affine_Fp3_Point(Right,Result,Result);
                                                                    csProjective:_Add_Projective_Fp3_Point(Right,Result,Result);
                                                                    csJacobian:_Add_Jacobian_Fp3_Point(Right,Result,Result);
                                                                  end;
                                                                 end;
                                         end;
     _Normalize_PF2_Point(Result,Result,CoordSys);
     if _IsNeg(left) then _Neg_Fp3_Point(Result,Result);
     end;
end;

       {**********   Multiply a scalar with an Fp Point *******}
procedure _Mul_NAF_Fp3_FpPoint(const Left:LInt;const Right:Fp3Point;var Result:Fp3Point;CoordSys:CoordinatesSystem=csAffine);
var i:integer;
    loop:LIntArrayForm;
    NegLeft:Fp3Point;
begin
Result.CurveParams:=Right.CurveParams;
_Neg_Fp3_Point(Right,NegLeft);
if (Right.Infinity)or (Left.Data.i16[-2]=0) then Result.Infinity:=true
else begin
     Result.Infinity:=true;
     if Left=2 then _Double_Jacobian_Fp3_Point(Right,Result)
     else begin
          Loop:=Left.ToNafArray;
          for i:=Length(Loop)-1 downto 0 do begin
                                            _Double_Fp3_Point(Result,Result,CoordSys);
                                            if Loop[i]=1 then
                                                _Add_Fp3_Point(Right,Result,Result,CoordSys);
                                            if Loop[i]=-1 then
                                                _Add_Fp3_Point(NegLeft,Result,Result,CoordSys);
                                            end;
          _Normalize_PF2_Point(Result,Result,CoordSys);
          end;
     end;
end;


       {**********   Compute Frobenius Map for an Fp3 Point *******}
procedure _Frobenius_Map_i(const Value:Fp3Point;pow:integer;var Result:Fp3Point);
var i:integer;
    tmp:Fp6Int;
    u:Fp3Int;
begin
if Value.Infinity then Result.Infinity:=true
else begin
     Result.SetCurveParams(Value.CurveParams);
     Result.X:=Value.X;
     Result.Y:=Value.Y;
     Result.Z:=Value.Z;
     if Value.CurveParams.Family=cfMNT then begin
                                            tmp.SetTowerParams(Value.CurveParams.TowerParam);
                                            u.SetFieldParams(Value.CurveParams.FieldParam);
                                            u.SetToZero;
                                            u.b:=1;

                                            /// Untwist
                                            Result.Y:=Lint(Value.CurveParams.FieldParam.Beta).InversModulo(Value.CurveParams.P).Sqr*Result.Y*u;
                                            tmp.a.a:=0;tmp.a.b:=Result.Y.b;
                                            tmp.b.a:=Result.Y.a;tmp.b.b:=0;
                                            tmp.c.a:=0;tmp.c.b:=Result.Y.c;
                                            ///*///
                                            for i:=0 to pow-1 do begin
                                                                 Result.X:=Result.X.toPowerP;
                                                                 _Pow_FP6_P(tmp,tmp);
                                                                 end;
                                            // Twist Back
                                            Result.Y.a:=tmp.b.a;Result.Y.b:=tmp.a.b;Result.Y.c:=tmp.c.b;
                                            Result.Y:=Lint(Value.CurveParams.FieldParam.Beta).Sqr*Result.Y*u.Inverse;
                                            ///*///
                                            end
     else begin
     for i:=0 to pow-1 do begin
                          _Mul_Fp3(Result.CurveParams.FrobeniusMapConstX_FP3,Result.X.toPowerP,Result.X);
                          _Mul_Fp3(Result.CurveParams.FrobeniusMapConstY_Fp3,Result.Y.toPowerP,Result.Y);
                          end;
          end;
    end;
end;

   {**********   Test if an Fp3 Point is on the Curve ***********}
function _Is_OnCurve_Fp3Point(value:Fp3Point):boolean;
var tmp,tmp1:Fp3Int;
begin
if value.Infinity then Result:=true
else begin
     _Sqr_Fp3(Value.X,tmp);
     _Add_Fp3(tmp,Value.CurveParams.AtwFp3,tmp);
     _Mul_Fp3(tmp,Value.X,tmp1);
     _Add_Fp3(tmp1,Value.CurveParams.AtwFp3,tmp);
     _Sqr_Fp3(Value.Y,tmp1);
    Result:=(tmp.IsASquare) and (tmp1=tmp);
    end;
end;

        {**********   Negation of an Fp3 Point *************}
procedure _Neg_Fp3_Point(Value:Fp3Point;var Result:Fp3Point);
begin
  Result:=Value;
  _Neg_Fp3(Result.Y,Result.Y);
end;

       {**********   Test if two Fp3 Points are equal *******}
function _Are_Equals_Fp3Points(left,right:Fp3Point):boolean;
begin
Result:=(_Compare_Fp3(Left.X,Right.X)=0)and(_Compare_Fp3(Left.Y,Right.Y)=0);
end;

      {**********   Evaluate an Fp3 Point from X value ******}
procedure _Evaluate_Fp3_Point(value:Fp3Int; var Result:Fp3Point);
var tmp,tmp1:Fp3Int;
begin
Result.X:=Value;
Result.Infinity:=false;
_Sqr_Fp3(Result.X,tmp);
_Add_Fp3(tmp,Result.CurveParams.AtwFp3,tmp);
_Mul_Fp3(tmp,Result.X,tmp1);
_Add_Fp3(tmp1,Result.CurveParams.BtwFp3,tmp);
if tmp.IsASquare then begin
                      _Sqrt_Fp3(tmp,Result.Y);
                      _Neg_Fp3(Result.Y,tmp1);
                      if _Compare_Fp3(Result.Y,tmp1)=1 then Result.Y:=tmp1;
                      Result.Z.a:=1;
                      Result.Z.b:=0;
                      Result.Z.c:=0;
                      end
else Result.Infinity:=true;
end;

{*******************************************************************************}
      ///  Definitions of an Fp3 Point operators and functions
{*******************************************************************************}

{*******************************************************************************}
procedure Fp3Point.SetCurveParams(Value: PtrCurveParams;ComputLmd:boolean=false);
begin
CurveParams:=Value;
if CurveParams.TowerParam3<>nil then LineAtP.SetTowerParams (Value.TowerParam3);
X.SetFieldParams(Value.FieldParam);
Y.SetFieldParams(Value.FieldParam);
Z.SetFieldParams(Value.FieldParam);
ComputeLigneValue:=ComputLmd;
ComputeLigneAtFp6:=ComputLmd;
Infinity:=false;
end;
{*******************************************************************************}
function Fp3Point.CompressToArray: TBytes;
var L:integer;
begin
if Infinity then begin
                 Setlength(result,0);
                 exit;
                 end;
L:=X.Field.p.Data.i32[-1]*4*3;
Setlength(Result,L);
Move(X.a.Data.i8[0],Result[0],Length(Result)div 3);
Move(X.b.Data.i8[0],Result[Length(Result)div 3],Length(Result)div 3 );
Move(X.c.Data.i8[0],Result[(Length(Result)div 3) shl 1],Length(Result)div 3);
end;
{*******************************************************************************}
procedure Fp3Point.DeCompressFromArray(a: TBytes);
begin
Move(a[0],X.a.Data.i8[0],Length(a) div 3);
Move(a[Length(a) div 3],X.b.Data.i8[0],Length(a) div 3);
Move(a[(Length(a) div 3)shl 1],X.c.Data.i8[0],Length(a) div 3);
X.a.Data.i16[-2]:=Length(a) shr 3;
X.a.Data.i16[-1]:=0;
X.b.Data.i16[-2]:=Length(a) shr 3;
X.b.Data.i16[-1]:=0;
X.c.Data.i16[-2]:=Length(a) shr 3;
X.c.Data.i16[-1]:=0;
SetPointValue(X);
end;
{*******************************************************************************}
function Fp3Point.FrobeniusMap_i(i:integer): Fp3Point;
begin
_Frobenius_Map_i(Self,i,Result);
end;
{*******************************************************************************}
procedure Fp3Point.InitFromStrings(Xa,Xb,Xc,Ya,Yb,Yc:String;Za:String='1';Zb:String='0';Zc:String='0');
begin
X.a:=Xa;
X.b:=Xb;
X.c:=Xc;
Y.a:=Ya;
Y.b:=Yb;
Y.c:=Yc;
Z.a:=Za;
Z.b:=Zb;
Z.c:=Zc;
Infinity:=false;
end;
{*******************************************************************************}
function Fp3Point.IsOnTheCurve: boolean;
begin
Result:=_Is_OnCurve_Fp3Point(Self)
end;
{*******************************************************************************}
procedure Fp3Point.SetAsRandomPoint;
begin
Infinity:=false;
Randomize;
repeat
  GetRandomLIntLowerThan(X.a,CurveParams.P);
  GetRandomLIntLowerThan(X.b,CurveParams.P);
  SetPointValue(X);
  until not Infinity;
end;
{*******************************************************************************}
procedure Fp3Point.SetAsTorsionFromHash(hash: TBytes);
var i,j:integer;
    h1,h2:TBytes;
    tmp:Fp3Point;
begin
j:=0;
repeat
  h1:=SHA256BytesHash(hash);
  h1[0]:=j;
  h2:=SHA256BytesHash(h1);
  Infinity:=false;
  X.a.Data.i32[length(h1) div 4]:=0;
  for i:=0 to length(h1)-1 do X.a.Data.i8[i]:=h1[i];
  X.a.Data.i16[-1]:=0;
  X.a.Data.i16[-2]:=(length(h1) div 4)+1;
  X.a:=X.a mod CurveParams.p;
  X.b.Data.i32[length(h1) div 4]:=0;
  for i:=0 to length(h1)-1 do X.b.Data.i8[i]:=h2[i];
  X.b.Data.i16[-1]:=0;
  X.b.Data.i16[-2]:=(length(h1) div 4)+1;
  X.b:=X.b mod CurveParams.p;
  X.c:=0;
  CurrentXp:=0;
  CurrentYp:=0;
  SetPointValue(X);
  if infinity then j:=j+1;
  until not Infinity;
tmp.SetCurveParams(Self.CurveParams);
//_Mul_Fp_Fp3Point(CurveParams.Htw,Self,tmp,csAffine);
if CurveParams.Family=cfKSS18 then _Fast_KSS18_Mul_Htw_Fp3Point(Self,tmp)
else if CurveParams.Family=cfMNT then  _Fast_MNT_Mul_Htw_Fp3Point(Self,tmp)
else _Mul_Fp_Fp3Point(CurveParams.Htw,Self,tmp,csAffine); // for MNT curves only affine coordinates are implemented
  Self:=tmp;
end;

{*******************************************************************************}
procedure Fp3Point.SetAsTorsionFromString(s: String);
begin
Self.SetAsTorsionFromHash(SHA256StringHash(s));
end;

{*******************************************************************************}
procedure Fp3Point.SetAsRandomTorsionPoint;
var tmp:Fp3Point;
  s:string;
begin
tmp.SetCurveParams(CurveParams);
Infinity:=false;
Randomize;
repeat
  GetRandomLIntLowerThan(X.a,CurveParams.P);
  GetRandomLIntLowerThan(X.b,CurveParams.P);
  GetRandomLIntLowerThan(X.c,CurveParams.P);
  SetPointValue(X);
  if not Infinity then begin
                       if CurveParams.Family=cfKSS18 then _Fast_KSS18_Mul_Htw_Fp3Point(Self,tmp)
                       else if CurveParams.Family=cfMNT then  _Fast_MNT_Mul_Htw_Fp3Point(Self,tmp)
                       else _Mul_Fp_Fp3Point(CurveParams.Htw,Self,tmp,csAffine); // for MNT curves only affine coordinates are implemented
                       Self:=tmp;
                       end;
  until (not Infinity) ;
  //_Mul_Fp_Fp3Point(CurveParams.Rtw,Self,tmp,csAffine);
//if not tmp.Infinity then showmessage('problem! G2');
end;
{*******************************************************************************}
procedure Fp3Point.SetPairingPointCoordinates(PtX, PtY: LInt);
begin
CurrentXp:=PtX;
CurrentYp:=PtY;
end;
{*******************************************************************************}
procedure Fp3Point.SetPointValue(Value: Fp3Int);
begin
_Evaluate_Fp3_Point(Value,Self);
end;

{*******************************************************************************}
procedure Fp3Point.SetToDefaultGenerator;
var tmp:Fp3Point;
begin
tmp.X:=CurveParams.TwistBasePointX_Fp3;
tmp.Y:=CurveParams.TwistBasePointy_Fp3;
tmp.Z.a:=1;
tmp.Z.b:=0;
tmp.Z.c:=0;
tmp.Infinity:=false;
tmp.SetCurveParams(CurveParams);
Tmp.ComputeLigneAtFp6:=false;
//Self:=CurveParams.Htw*tmp;
if CurveParams.Family=cfKSS18 then _Fast_KSS18_Mul_Htw_Fp3Point(tmp,Self)
else if CurveParams.Family=cfMNT then  _Fast_MNT_Mul_Htw_Fp3Point(tmp,Self)
else _Mul_Fp_Fp3Point(CurveParams.Htw,tmp,self,csAffine); // for MNT curves only affine coordinates are implemented

end;

{*******************************************************************************}
function Fp3Point.ToHexString: String;
begin
if Infinity then Result:='Infinity' else Result:='('+X.ToHexString+','+Y.ToHexString+')';
end;
{*******************************************************************************}
function Fp3Point.ToDecimalString: String;
begin
if Infinity then Result:='Infinity' else Result:='('+X.ToDecimalString+','+Y.ToDecimalString+')';
end;
{*******************************************************************************}
class operator Fp3Point.Add(Left, Right: Fp3Point): Fp3Point;
begin
_Add_Affine_Fp3_Point(Left,Right,Result);
end;
{*******************************************************************************}
class operator Fp3Point.Multiply(Left: LInt; Right{Q}: Fp3Point): Fp3Point;
begin
_Mul_Fp_Fp3Point(Left,Right,Result);
end;
{*******************************************************************************}
class operator Fp3Point.Subtract(Left, Right: Fp3Point): Fp3Point;
begin
_Sub_Fp3_Point(left,right,result);
end;
{*******************************************************************************}
class operator Fp3Point.Negative(const Value: Fp3Point): Fp3Point;
begin
_Neg_Fp3_Point(Value,Result);
end;
{*******************************************************************************}
class operator Fp3Point.Equal(const Left, Right: Fp3Point): Boolean;
begin
Result:=_Are_Equals_Fp3Points(Left,Right);
end;
{*******************************************************************************}
class operator Fp3Point.NotEqual(const Left, Right: Fp3Point): Boolean;
begin
Result:=not _Are_Equals_Fp3Points(Left,Right);
end;

end.
