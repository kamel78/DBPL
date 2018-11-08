unit Fp6Arithmetic;

interface
uses System.SysUtils,VCL.dialogs, LargeIntegers, Fp2Arithmetic,Fp3Arithmetic;


{*******************************************************************************************
   Développer par Faraoun Kamel Mohamed
   Université Djilali Liabes -Sidi Bel Abbes - Algérie
   kamel_mh@yahoo.fr

   Arithmetic computation over Fp6. Fp6 is the sextic extention of Fp (Tower extention of order 3 for FP2)
   with respect to the irriductible polynômial V^3-Sigma (V^3=Sigma). Elements are in the form a+b*V+c*V^2
   where a,b and c are from Fp2
********************************************************************************************}

type

  PtrTowerParams12=^TowerParams12;
  TowerParams12=record
              FieldParam:PtrFieldParams;    // Parmetres of the Base Field
              Sigma:Fp2Int;                 // Parametre of the Cubic Extention
              inv_2:Fp2Int;
              pmod8:Lint;                   //used for BLS Curves
              FrobeniusP_Const,FrobeniusP2_Const,
              FrobeniusP3_Const:array[0..4]of Fp2Int;  // Frobenius Constants for PreComputation of f^p, f^p^2, f^p^3 (BN Curves)
              FrobeniusPi3xSigmaSqr:Fp2Int;       // Frobenius Constants for PreComputation of Fp4_Power_P (BLS24 Curves)
              end;

  Fp6Int=record
  a,b,c:Fp2Int;
  Tower:PtrTowerParams12;
  public
    class operator Add(const Left, Right: Fp6Int): Fp6Int;
    class operator Subtract(const Left, Right: Fp6Int): Fp6Int;
    class operator Multiply(const Left, Right: Fp6Int): Fp6Int;
    class operator Multiply(const Left:Fp2Int; Right: Fp6Int): Fp6Int;
    class operator Negative(const Value: Fp6Int): Fp6Int;
    class operator Equal(const Left, Right: Fp6Int): Boolean;
    class operator NotEqual(const Left, Right: Fp6Int): Boolean;
    function ToHexString: string;
    function ToDecimalString: string;
    function ToByteArray:Tbytes;
    function Inverse:Fp6Int;
    function Pow(Exponent: LInt): Fp6Int;
    function IsZero:boolean;
    function IsOne:boolean;
    function Sqr:Fp6Int;
    function Sqrt:Fp6Int;
    function MultiplyByGamma:Fp6Int;
    procedure SetToRandom;
    procedure SetToOne;
    procedure SetToZero;
    procedure SetTowerParams(TParam:PtrTowerParams12);
    procedure SetFromStrings(a,b,c:String);
    function toPowerP: Fp6Int;
    function IsASquare: Boolean;
  end;

  procedure _Add_FP6(const Left,Right:Fp6Int; var Result:Fp6Int);
  procedure _Sub_FP6(Const Left,Right:Fp6Int;var Result: Fp6Int);
  procedure _Mul_FP6(const Left,Right:Fp6Int;var Result:Fp6Int);
  procedure _Mul_FP2_FP6(const Left: Fp2Int; Right: Fp6Int; var Result: Fp6Int);
  procedure _Mul_FP_FP6(const Left: LInt; Right: Fp6Int; var Result: Fp6Int);
  procedure _Mul_FP6_By_V2(const value: Fp6Int; var Result: Fp6Int);
  procedure _Sqr_FP6(const value: Fp6Int; var Result: Fp6Int);
  function _Equals_FP6(const Left,Right:Fp6Int):boolean;
  procedure _Inv_FP6(const value:Fp6Int;var Result:Fp6Int);
  procedure _Neg_FP6(const Value: Fp6Int; var Result: Fp6Int);
  procedure _Pow_FP6(const value:Fp6Int;Exponent:LInt;var Result:Fp6Int);
  function _FP6_To_DecimalString(const value:Fp6Int):string;
  function _FP6_To_HexString(const value:Fp6Int):string;
  procedure _FP6_From_Strings(const valueA,valueB,valueC:String;var Result:Fp6Int);
  procedure _Sparse1_Mul_FP6(const Left,Right:Fp6Int;var Result:Fp6Int);
  procedure _Sparse2_Mul_FP6(const Left,Right:Fp6Int;var Result:Fp6Int);
  procedure _Sparse3_Mul_FP6(const Left,Right:Fp6Int;var Result:Fp6Int);
  procedure _Sparse_Add_FP6(const Left,Right:Fp6Int; var Result:Fp6Int);
  procedure _Nomod_Sub_FP6(Const Left,Right:Fp6Int;var Result: Fp6Int);
  procedure _Pow_FP6_P(const value:Fp6Int;var Result:Fp6Int);
  procedure _Pow_FP6_P_i(const value:Fp6Int; power: integer;var Result:Fp6Int);
  function _Is_Fp6_ASquare(const value:Fp6Int):boolean;
  procedure _Sqrt_Fp6(const Value:Fp6Int;var Result:Fp6Int);

implementation

{*******************************************************************************}
            ///      Procedures for FP6 Arithmetic
{*******************************************************************************}

    {**********   Add two Sparse FP6 integers *****************}
    //  Left is in Fp and Right.C is 0 (for D-type Twist only !)
procedure _Sparse_Add_FP6(const Left,Right:Fp6Int; var Result:Fp6Int);
begin
Result.Tower:=Left.Tower;
_Add_FP2(Left.a,Right.a,Result.a);
Result.b:=Right.b;
Result.c.a.Data.i16[-2]:=0;
Result.c.b.Data.i16[-2]:=0;
end;

    {**********   Add two FP6 integers *****************}
procedure _Add_FP6(const Left,Right:Fp6Int; var Result:Fp6Int);
begin
Result.Tower:=Left.Tower;
_Add_FP2(Left.a,Right.a,Result.a);
_Add_FP2(Left.b,Right.b,Result.b);
_Add_FP2(Left.c,Right.c,Result.c);
end;

    {**********   Substract two FP6 integers ***********}
procedure _Sub_FP6(Const Left,Right:Fp6Int;var Result: Fp6Int);
begin
Result.Tower:=Left.Tower;
_Sub_FP2(Left.a,Right.a,Result.a);
_Sub_FP2(Left.b,Right.b,Result.b);
_Sub_FP2(Left.c,Right.c,Result.c);
end;

    {**********   Substract two FP6 integers ***********}
     // Without Modular reduction (for optimization in some cases)
procedure _Nomod_Sub_FP6(Const Left,Right:Fp6Int;var Result: Fp6Int);
begin
Result.Tower:=Left.Tower;
_Nomod_Sub_FP2(Left.a,Right.a,Result.a);
_Nomod_Sub_FP2(Left.b,Right.b,Result.b);
_Nomod_Sub_FP2(Left.c,Right.c,Result.c);
end;

    {********** Sparse Multiply two FP6 integers *****************}
    // Sparse_1 Multiplication Using Karatsuba
    // Right is in FP2 (Alwyas Sparse)  , for D-type Twist only !
procedure _Sparse1_Mul_FP6(const Left,Right:Fp6Int;var Result:Fp6Int);
begin
Result.Tower:=Left.Tower;
_Mul_FP2(Left.a,Right.a,Result.a);
_Mul_FP2(Left.b,Right.a,Result.b);
_Mul_FP2(Left.c,Right.a,Result.c);
end;

    {********** Sparse Multiply two FP6 integers *****************}
    //Sparse 1 Multiplication Using Karatsuba Multiplication
    // Right.C=0 and right.a=0 and right.b.b=0 (only right.b.a<>0). for M-type Twist only !
procedure _Sparse3_Mul_FP6(const Left,Right:Fp6Int;var Result:Fp6Int);
var tmp:array[0..1] of Fp2Int;
begin
Result.SetTowerParams(Left.Tower);
_Add_FP2(Left.b,Left.c,tmp[0]);
_Add_FP2(Left.a,Left.b,tmp[1]);
_Mul_FP2(Left.b,Right.b,Result.c);
_Mul_FP2(tmp[0],right.b,Result.a);
_Nomod_Sub_FP2(Result.a,Result.c,Result.a);
_Mul_FP2(Result.a,Left.Tower^.sigma,Result.a);
_Mul_FP2(tmp[1],right.b,Result.b);
_Sub_FP2(Result.b,Result.c,Result.b);
end;

        {********** Sparse Multiply two FP6 integers *****************}
        // Sparse-2 Multiplication Using Karatsuba Multiplication
        /// Right.c is always 0 (Partially Sparse)   , for D-type Twist only !)
procedure _Sparse2_Mul_FP6(const Left,Right:Fp6Int;var Result:Fp6Int);
var tmp:array[0..6] of Fp2Int;
begin
Result.SetTowerParams(Left.Tower);
_Mul_FP2(Left.a,Right.a,tmp[0]);
_Mul_FP2(Left.b,Right.b,tmp[1]);
_Nomod_Add_FP2(Left.b,Left.c,tmp[3]);
_NoMod_Add_FP2(Left.a,Left.b,tmp[4]);
_NoMod_Add_FP2(Right.a,Right.b,tmp[5]);
_Nomod_Add_FP2(Left.a,Left.c,tmp[6]);
_Mul_FP2(tmp[3],right.b,Result.a);
_Nomod_Sub_FP2(Result.a,tmp[1],Result.a);
_Mul_FP2(Result.a,Left.Tower^.sigma,Result.a);
_Add_FP2(tmp[0],Result.a,Result.a);
_Mul_FP2(tmp[4],tmp[5],Result.b);
_Nomod_Sub_FP2(Result.b,tmp[0],Result.b);
_Sub_FP2(Result.b,tmp[1],Result.b);
_Mul_FP2(tmp[6],right.a,Result.c);
_Nomod_Sub_FP2(Result.c,tmp[0],Result.c);
_Add_FP2(Result.c,tmp[1],Result.c);
end;

        {********** Multiply two FP6 integers *****************}
procedure _Mul_FP6(const Left,Right:Fp6Int;var Result:Fp6Int);
var tmp:array[0..9] of Fp2Int;
begin
        { Multiplication Using Karatsuba Multiplication  }
Result.SetTowerParams(Left.Tower);
_Mul_FP2(Left.a,Right.a,tmp[0]);
_Mul_FP2(Left.b,Right.b,tmp[1]);
_Mul_FP2(Left.c,Right.c,tmp[2]);
_Add_FP2(Left.b,Left.c,tmp[3]);
_Add_FP2(Right.b,Right.c,tmp[4]);
_Add_FP2(Left.a,Left.b,tmp[5]);
_Add_FP2(Right.a,Right.b,tmp[6]);
_Add_FP2(Left.a,Left.c,tmp[7]);
_Add_FP2(Right.a,Right.c,tmp[8]);
_Mul_FP2(tmp[3],tmp[4],Result.a);
_Nomod_Sub_FP2(Result.a,tmp[1],Result.a);
_Sub_FP2(Result.a,tmp[2],Result.a);
_Mul_FP2(Result.a,Left.Tower^.sigma,Result.a);
_Add_FP2(tmp[0],Result.a,Result.a);
_Mul_FP2(tmp[5],tmp[6],Result.b);
_Nomod_Sub_FP2(Result.b,tmp[0],Result.b);
_Sub_FP2(Result.b,tmp[1],Result.b);
_Mul_FP2(tmp[2],Left.Tower^.sigma,tmp[9]);
_Add_FP2(Result.b,tmp[9],Result.b);
_Mul_FP2(tmp[7],tmp[8],Result.c);
_Nomod_Sub_FP2(Result.c,tmp[0],Result.c);
_Nomod_Sub_FP2(Result.c,tmp[2],Result.c);
_Add_FP2(Result.c,tmp[1],Result.c);
end;

        {********** Multiply FP2 with FP6 *****************}
procedure _Mul_FP2_FP6(const Left: Fp2Int; Right: Fp6Int; var Result: Fp6Int);
begin
Result.SetTowerParams(Right.Tower);
_Mul_FP2(Right.a,Left,Result.a);
_Mul_FP2(Right.b,Left,Result.b);
_Mul_FP2(Right.c,Left,Result.c);
end;

        {********** Multiply FP with FP6 *****************}
procedure _Mul_FP_FP6(const Left: LInt; Right: Fp6Int; var Result: Fp6Int);
begin
Result.SetTowerParams(Right.Tower);
_Mul_Fp_FP2(Left,Right.a,Result.a);
_Mul_FP_FP2(Left,Right.b,Result.b);
_Mul_Fp_FP2(Left,Right.c,Result.c);
end;

        {********** Multiply FP2 with FP6 *****************}
procedure _Mul_FP6_By_V2(const value: Fp6Int; var Result: Fp6Int);  // Value should be different than Result
var tmp:Fp2Int;
begin
Result.SetTowerParams(Value.Tower);
_Mul_FP2(Value.c,Value.Tower^.Sigma,tmp);
Result.c:=Value.b;
Result.b:=Value.a;
Result.a:=tmp;
end;
        {********** Get Square of an FP6 *****************}
procedure _Sqr_FP6(const value: Fp6Int; var Result: Fp6Int);
var tmp: array [0..4] of Fp2Int;
begin
Result.SetTowerParams(Value.Tower);
      {       Using Karatsuba Squering    }
_Sqr_FP2(Value.a,tmp[0]);
_Sqr_FP2(Value.b,tmp[1]);
_Sqr_FP2(Value.c,tmp[2]);
_Add_FP2(Value.a,Value.b,tmp[3]);
_Add_FP2(Value.a,Value.c,tmp[4]);
_Add_FP2(Value.b,Value.c,Result.a);
_Sqr_FP2(Result.a,Result.a);
_Sub_FP2(Result.a,tmp[1],Result.a);
_Sub_FP2(Result.a,tmp[2],Result.a);
_Mul_FP2(value.Tower^.Sigma,Result.a,Result.a);
_Add_FP2(Result.a,tmp[0],Result.a);
_Sqr_FP2(tmp[3],tmp[3]);
_Sub_FP2(tmp[3],tmp[0],tmp[3]);
_Sub_FP2(tmp[3],tmp[1],tmp[3]);
_Mul_FP2(Value.Tower^.Sigma,tmp[2],Result.b);
_Add_FP2(Result.b,tmp[3],Result.b);
_Sqr_FP2(tmp[4],Result.c);
_Sub_FP2(Result.c,tmp[0],Result.c);
_Sub_FP2(Result.c,tmp[2],Result.c);
_Add_FP2(Result.c,tmp[1],Result.c);
end;

        {**********   Compare two FP6 integers *************}
function _Equals_FP6(const Left,Right:Fp6Int):boolean;
begin
result:=(_Compare_FP2(Left.a,Right.a)=0) and (_Compare_FP2(Left.b,Right.b)=0)and (_Compare_FP2(Left.c,Right.c)=0);
end;

        {********** Get inverse of an FP6 *****************}
procedure _Inv_FP6(const value:Fp6Int;var Result:Fp6Int);
var t:array[0..9] of Fp2Int;
begin
Result.Tower:=Value.Tower;
_Sqr_FP2(value.a,t[0]);
_Sqr_FP2(value.b,t[1]);
_Sqr_FP2(value.c,t[2]);
_Mul_FP2(value.a,value.b,t[3]);
_Mul_FP2(value.a,value.c,t[4]);
_Mul_FP2(value.b,value.c,t[5]);
_Mul_FP2(value.Tower^.Sigma,t[5],t[7]);
_Sub_FP2(t[0],t[7],t[7]);
_Mul_FP2(value.Tower^.Sigma,t[2],t[8]);
_Sub_FP2(t[8],t[3],t[8]);
_Sub_FP2(t[1],t[4],t[9]);
_Mul_FP2(value.a,t[7],t[6]);
_Mul_FP2(value.Tower^.Sigma,value.c,t[0]);
_Mul_FP2(t[0],t[8],t[0]);
_Add_FP2(t[6],t[0],t[6]);
_Mul_FP2(value.Tower^.Sigma,value.b,t[0]);
_Mul_FP2(t[0],t[9],t[0]);
_Add_FP2(t[6],t[0],t[6]);
_Inv_FP2(t[6],t[6]);
_Mul_FP2(t[7],t[6],Result.a);
_Mul_FP2(t[8],t[6],Result.b);
_Mul_FP2(t[9],t[6],Result.c);
end;
       {****** Test if an Fp6 is a Sequare *******}
function _Is_Fp6_ASquare(const value:Fp6Int):boolean;
var t1,t2:fp3Int;
begin
if value.IsZero then result:=true
else begin
     if value.b.a.IsZero and value.a.b.IsZero and value.c.b.IsZero then result:=true
else begin
     t1.SetFieldParams(value.Tower.FieldParam);
     t2.SetFieldParams(value.Tower.FieldParam);
     _HCopy_LInt(value.a.a,t1.a);
     _HCopy_LInt(value.c.a,t1.b);
     _HCopy_LInt(value.b.b,t1.c);
     _HCopy_LInt(value.b.a,t2.a);
     _HCopy_LInt(value.a.b,t2.b);
     _HCopy_LInt(value.c.b,t2.c);
     _Sqr_Fp3(t1,t1);
     _Sqr_Fp3(t2,t2);
     _Mul_Fp3_u(t2,t2);
     _Sub_Fp3(t1,t2,t1);
     result:=t1.IsASquare;
     end;
     end;
end;

        {********** Get Square root of an Fp6 ************}
        /// Convert from 2-3 to 3-2, find the root then return to 2-3
procedure _Sqrt_Fp6(const Value:Fp6Int;var Result:Fp6Int);
var t1,t2,s,a:Fp3Int;
begin
Result.SetTowerParams(Value.Tower);
if not _Is_Fp6_ASquare(Value) then  raise ERangeError.Create( 'This FP6 element is not a Square .....')
else begin
     if value.IsZero then result.SetToZero
     else begin
          t1.SetFieldParams(value.Tower.FieldParam);
          t2.SetFieldParams(value.Tower.FieldParam);
          a.SetFieldParams(value.Tower.FieldParam);
          s.SetFieldParams(value.Tower.FieldParam);
          // Convert from towering 2-3 to towering 3-2
          _HCopy_LInt(value.a.a,t1.a);
          _HCopy_LInt(value.c.a,t1.b);
          _HCopy_LInt(value.b.b,t1.c);
          _HCopy_LInt(value.b.a,t2.a);
          _HCopy_LInt(value.a.b,t2.b);
          _HCopy_LInt(value.c.b,t2.c);
          //
          if t2.IsZero then begin
                            if t1.IsASquare then begin
                                                 _Sqrt_Fp3(t1,t1);
                                                 // Convert from towering 3-2 to towering 2-3
                                                 _HCopy_LInt(t1.a,Result.a.a);
                                                 _HCopy_LInt(t1.b,Result.c.a);
                                                 _HCopy_LInt(t1.c,Result.b.b);
                                                 Result.b.a:=0;
                                                 Result.a.b:=0;
                                                 Result.c.b:=0;;
                                                 end
                            else begin
                                 _Div_Fp3_u(t1,t1);
                                 _Sqr_Fp3(t1,t1);
                                 // Convert from towering 3-2 to towering 2-3
                                 Result.a.a:=0;
                                 Result.c.a:=0;
                                 Result.b.b:=0;
                                 _HCopy_LInt(t1.a,Result.b.a);
                                 _HCopy_LInt(t1.b,Result.a.b);
                                 _HCopy_LInt(t1.c,Result.c.b);
                                 end;
                            end
          else begin
               _Sqr_Fp3(t1,a);
               _Sqr_Fp3(t2,s);
               _Mul_Fp3_u(s,s);
               _Sub_Fp3(a,s,a);
               _Sqrt_Fp3(a,s);
               if s.IsZero then Result.SetToZero
               else begin
                    _Add_Fp3(t1,s,a);
                    _Mul_FP_Fp3(Value.Tower.FieldParam.inv_2,a,a);
                    if a.IsASquare then _Sqrt_Fp3(a,a)
                    else begin
                         _Sub_Fp3(t1,s,a);
                         _Mul_FP_Fp3(Value.Tower.FieldParam.inv_2,a,a);
                         _Sqrt_Fp3(a,a);
                         if a.IsZero then begin
                                          Result.SetToZero;
                                          exit;
                                          end;
                         end;
                         // Convert from towering 3-2 to towering 2-3
                         _HCopy_LInt(a.a,Result.a.a);
                         _HCopy_LInt(a.b,Result.c.a);
                         _HCopy_LInt(a.c,Result.b.b);
                         _Mul_FP_Fp3(Value.Tower.FieldParam.inv_2,t2,s);
                         _Mul_Fp3(a.Inverse,s,s);
                         _HCopy_LInt(s.a, Result.b.a);
                         _HCopy_LInt(s.b,Result.a.b);
                         _HCopy_LInt(s.c,Result.c.b);
                    end;
               end;
          end
     end;
end;

        {********** Get negative of an FP6 *****************}
procedure _Neg_FP6(const Value: Fp6Int; var Result: Fp6Int);
begin
Result.Tower:=Value.Tower;
_Neg_FP2(Value.a,Result.a);
_Neg_FP2(Value.b,Result.b);
_Neg_FP2(Value.c,Result.c);
end;

        {********** Raise an FP6 to a LInt power ************}
procedure _Pow_FP6(const value:Fp6Int;Exponent:LInt;var Result:Fp6Int);
var i:word;
    tmp:Fp6Int;
begin
Result.Tower:= Value.Tower;
Result.a:=Value.a;
Result.b:=Value.b;
Result.c:=Value.c;
tmp:=Value;
if (Exponent.Absolute>1) then
for i:=Exponent.BitLength-2 downto 0 do begin
                                        _Sqr_FP6(Result,Result);
                                        if _Is_BitSet_At(Exponent,i) then _Mul_FP6(Result,tmp,Result);
                                        end;
end;


        {********** Raise an Fp6 to " p " power ************}
        /// use Precomputed Frobenius Constants (For use in MNT/KSS36 Curves !)
procedure _Pow_FP6_P(const value:Fp6Int;var Result:Fp6Int);
begin
Result.Tower:=Value.Tower;
Result.a.SetFieldParams(Value.Tower.FieldParam);
Result.b.SetFieldParams(Value.Tower.FieldParam);
Result.c.SetFieldParams(Value.Tower.FieldParam);
Result.a:=Value.a.Conjugate;
Result.b:=Value.b.Conjugate;
Result.c:=Value.c.Conjugate;
_Mul_FP2(Result.b,Value.Tower^.FrobeniusP_Const[0],Result.b);
_Mul_FP2(Result.c,Value.Tower^.FrobeniusP_Const[1],Result.c);
end;

        {********** Raise an Fp6 to " p^i " power ************}
        /// use Precomputed Frobenius Constant  (For use in MNT/KSS36 Curves exclusivly !)
procedure _Pow_FP6_P_i(const value:Fp6Int; power: integer;var Result:Fp6Int);
var i:integer;
begin
Result.Tower:=Value.Tower;
Result.a.SetFieldParams(Value.Tower.FieldParam);
Result.b.SetFieldParams(Value.Tower.FieldParam);
Result.c.SetFieldParams(Value.Tower.FieldParam);
Result := Value;
for i := 0 to power - 1 do _Pow_FP6_P(Result, Result);
end;


        {****** Convert an FP6 to a Decimal String **********}
function _FP6_To_DecimalString(const value:Fp6Int):string;
begin
if Value.c.IsZero then begin
                 if Value.b.IsZero then begin
                                  if Value.a.IsZero then Result:='0'
                                  else Result:=Value.a.ToDecimalString;
                                  end
                 else begin
                      if Value.b.IsOne then begin
                                      if Value.a.IsZero then Result:='v'
                                      else Result:='v +'+Value.a.ToDecimalString;
                                      end
                       else begin
                            if Value.a.IsZero then Result:='('+Value.b.ToDecimalString+') * v'
                            else Result:='('+Value.b.ToDecimalString+') * v +'+Value.a.ToDecimalString;
                            end;
                      end
                 end
else begin
     if Value.c.IsOne then begin
                     if Value.b.IsZero then begin
                                      if Value.a.IsZero then Result:='v^2'
                                      else Result:='v^2 + '+Value.a.ToDecimalString;
                                      end
                     else begin
                          if Value.b.IsOne then begin
                                          if Value.a.IsZero then Result:='v^2 + v'
                                          else Result:='v^2 + v +'+Value.a.ToDecimalString;
                                          end
                          else begin
                               if Value.a.IsZero then Result:='v^2 +('+Value.b.ToDecimalString+') * v'
                               else Result:='v^2 + ('+Value.b.ToDecimalString+') * v +'+Value.a.ToDecimalString;
                               end;
                         end
                    end
    else begin
         if Value.b.IsZero then begin
                          if Value.a.IsZero then Result:='('+Value.c.ToDecimalString+') * v^2'
                                      else Result:='('+Value.c.ToDecimalString+') * v^2 + '+Value.a.ToDecimalString;
                                      end
                          else begin
                               if Value.b.IsOne then begin
                                               if Value.a.IsZero then Result:='('+Value.c.ToDecimalString+') * v^2 +'+ 'v'
                                               else Result:='('+Value.c.ToDecimalString+') * v^2 +'+' v +'+Value.a.ToDecimalString;
                                               end
                               else begin
                                    if Value.a.IsZero then Result:='('+Value.c.ToDecimalString+') * v^2 +'+'('+Value.b.ToDecimalString+') * v'
                                    else Result:='('+Value.c.ToDecimalString+') * v^2 +'+'('+Value.b.ToDecimalString+') * v +'+Value.a.ToDecimalString;
                                    end;
                              end
                          end;
        end;
end;

        {****** Convert an FP6 to a Hexadecimal String *******}
function _FP6_To_HexString(const value:Fp6Int):string;
begin
if Value.c.IsZero then begin
                 if Value.b.IsZero then begin
                                  if Value.a.IsZero then Result:='0'
                                  else Result:=Value.a.ToHexString;
                                  end
                 else begin
                      if Value.b.IsOne then begin
                                      if Value.a.IsZero then Result:='v'
                                      else Result:='v +'+Value.a.ToHexString;
                                      end
                       else begin
                            if Value.a.IsZero then Result:='('+Value.b.ToHexString+') * v'
                            else Result:='('+Value.b.ToHexString+') * v +'+Value.a.ToHexString;
                            end;
                      end
                 end
else begin
     if Value.c.IsOne then begin
                     if Value.b.IsZero then begin
                                      if Value.a.IsZero then Result:='v^2'
                                      else Result:='v^2 + '+Value.a.ToHexString;
                                      end
                     else begin
                          if Value.b.IsOne then begin
                                          if Value.a.IsZero then Result:='v^2 + v'
                                          else Result:='v^2 + v +'+Value.a.ToHexString;
                                          end
                          else begin
                               if Value.a.IsZero then Result:='v^2 +('+Value.b.ToHexString+') * v'
                               else Result:='v^2 + ('+Value.b.ToHexString+') * v +'+Value.a.ToHexString;
                               end;
                         end
                    end
    else begin
         if Value.b.IsZero then begin
                          if Value.a.IsZero then Result:='('+Value.c.ToHexString+') * v^2'
                                      else Result:='('+Value.c.ToHexString+') * v^2 + '+Value.a.ToHexString;
                                      end
                          else begin
                               if Value.b.IsOne then begin
                                               if Value.a.IsZero then Result:='('+Value.c.ToHexString+') * v^2 +'+ 'v'
                                               else Result:='('+Value.c.ToHexString+') * v^2 +'+' v +'+Value.a.ToHexString;
                                               end
                               else begin
                                    if Value.a.IsZero then Result:='('+Value.c.ToHexString+') * v^2 +'+'('+Value.b.ToHexString+') * v'
                                    else Result:='('+Value.c.ToHexString+') * v^2 +'+'('+Value.b.ToHexString+') * v +'+Value.a.ToHexString;
                                    end;
                              end
                          end;
        end;
end;

        {****** Convert String to an FP2 (Decimal/Hex)*******}
procedure _FP6_From_Strings(const valueA,valueB,valueC:String;var Result:Fp6Int);
var i:integer;
    s1,s2,Val:string;
    valid:boolean;
begin
Result.a.SetFormString(valueA);
Result.b.SetFormString(valueB);
Result.c.SetFormString(valueC);
end;

{*******************************************************************************}
      ///  Definitions of an FP6 integer operators and functions
{*******************************************************************************}
class operator Fp6Int.Add(const Left, Right: Fp6Int): Fp6Int;
begin
_Add_FP6(Left,Right,Result);
end;

class operator Fp6Int.Subtract(const Left, Right: Fp6Int): Fp6Int;
begin
_Sub_FP6(Left,Right,Result);
end;

class operator Fp6Int.Multiply(const Left, Right: Fp6Int): Fp6Int;
begin
_Mul_FP6(Left,Right,Result);
end;

class operator Fp6Int.Multiply(const Left: Fp2Int;  Right: Fp6Int): Fp6Int;
begin
_Mul_FP2_FP6(Left,Right,Result);
end;

function Fp6Int.MultiplyByGamma: Fp6Int;
begin
_Mul_FP6_By_V2(Self,Result);
end;

function Fp6Int.Sqr: Fp6Int;
begin
_Sqr_FP6(Self,Result);
end;

function Fp6Int.Sqrt: Fp6Int;
begin
_Sqrt_Fp6(Self,Result);
end;

class operator Fp6Int.Equal(const Left, Right: Fp6Int): Boolean;
begin
result:=_Equals_FP6(Left,Right);
end;

function Fp6Int.Inverse: Fp6Int;
begin
_Inv_FP6(Self,Result);
end;

function Fp6Int.IsOne: boolean;
begin
Result:=(a.IsOne)and(b.IsZero)and(c.IsZero);
end;

function Fp6Int.IsZero: boolean;
begin
Result:=(a.IsZero)and(b.IsZero)and(c.IsZero);
end;

class operator Fp6Int.Negative(const Value: Fp6Int): Fp6Int;
begin
_Neg_FP6(Value,Result);
end;

class operator Fp6Int.NotEqual(const Left, Right: Fp6Int): Boolean;
begin
Result:=not _Equals_FP6(Left,Right)
end;

function Fp6Int.Pow(Exponent: LInt): Fp6Int;
begin
_Pow_FP6(Self,Exponent,Result);
end;

procedure Fp6Int.SetTowerParams(TParam:PtrTowerParams12);
begin
Tower:=TParam;
a.SetFieldParams(Tparam.FieldParam);
b.SetFieldParams(Tparam.FieldParam);
c.SetFieldParams(Tparam.FieldParam);
end;

procedure Fp6Int.SetToZero;
begin
a.SetToZero;
b.SetToZero;
c.SetToZero;
end;

procedure Fp6Int.SetFromStrings(a, b, c: String);
begin
_FP6_From_Strings(a,b,c,Self);
end;

procedure Fp6Int.SetToOne;
begin
a.SetToZero;
b.SetToZero;
c.SetToZero;
a.a:=1;
end;

procedure Fp6Int.SetToRandom;
begin
a.SetToRandom;
b.SetToRandom;
c.SetToRandom;
end;

function Fp6Int.ToByteArray: Tbytes;
var size:integer;
begin
size:=(a.a.BitLength mod 8) * 8;
Setlength(Result,size*6);
Move(a.a.Data.i8[0],Result[0],size);
Move(a.b.Data.i8[0],Result[size],size);
Move(b.a.Data.i8[0],Result[2*size],size);
Move(b.b.Data.i8[0],Result[3*size],size);
Move(c.a.Data.i8[0],Result[4*size],size);
Move(c.b.Data.i8[0],Result[5*size],size);
end;

function Fp6Int.toDecimalString: string;
begin
Result:=_FP6_To_DecimalString(Self)
end;

function Fp6Int.toHexString: string;
begin
Result:=_FP6_To_HexString(Self)
end;

function Fp6Int.toPowerP: Fp6Int;
begin
_Pow_FP6_P(Self,Result);
end;

function Fp6Int.IsASquare: Boolean;
begin
Result := _Is_FP6_ASquare(Self);
end;


end.
