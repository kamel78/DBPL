unit HashFunctions;

interface

uses IdGlobal,IdHashSha,IdSSLOpenSSL,System.SysUtils,System.Classes;


Type

  Digest256Bit=Record
             Data:array [0..7] of UInt32;
             function ToString:AnsiString;
             function ToByteArray:TBytes;
             End;

  TSHA256=class
    private
    Fm_buffer_size, Fm_block_size, Fm_hash_size: Int32;
    Fm_pos: Int32;
    Fm_buffer: Tbytes;
    Fm_processed_bytes: UInt64;
    Fm_state:Digest256Bit;
    function Feed(a_data: PByte; a_length_a_data: Int32;var a_start_index, a_length: Int32; var a_processed_bytes: UInt64): Boolean;
    procedure Finish;
    procedure Initialize;
    procedure TransformBlock(a_data: PByte; a_data_length: Int32;a_index: Int32);
    procedure TransformBytes(a_data: TBytes;a_index, a_length: Int32);
    function ComputeBytes(a_data: Tbytes): Digest256Bit;
    public
    constructor create;
    function HashString(S:String):Digest256Bit;
    function HashByteArray(B:Tbytes):Digest256Bit;
    function HashFile(filename:String):Digest256Bit;
  end;



function SHA256FileHash(_filename: string): TBytes;
function SHA256StringHash(_string: string): Tbytes;
function SHA256BytesHash(_bytes: Tbytes): Tbytes;

implementation

var Sha256Hasher:TSHA256;


function RotateRight32(a_value: UInt32; a_n: Int32): UInt32;
begin
a_n := a_n and 31;
result := (a_value shr a_n) or (a_value shl (32 - a_n));
end;

function ReverseBytesUInt64(value: UInt64): UInt64;
begin
result := (value and UInt64($00000000000000FF)) shl 56 or
    (value and UInt64($000000000000FF00)) shl 40 or
    (value and UInt64($0000000000FF0000)) shl 24 or
    (value and UInt64($00000000FF000000)) shl 8 or
    (value and UInt64($000000FF00000000)) shr 8 or
    (value and UInt64($0000FF0000000000)) shr 24 or
    (value and UInt64($00FF000000000000)) shr 40 or
    (value and UInt64($FF00000000000000)) shr 56;
end;

procedure ReadUInt64AsBytesLE(a_in: UInt64;a_out: Tbytes; a_index: Int32);
begin
  a_out[a_index] := Byte(a_in);
  a_out[a_index + 1] := Byte(a_in shr 8);
  a_out[a_index + 2] := Byte(a_in shr 16);
  a_out[a_index + 3] := Byte(a_in shr 24);
  a_out[a_index + 4] := Byte(a_in shr 32);
  a_out[a_index + 5] := Byte(a_in shr 40);
  a_out[a_index + 6] := Byte(a_in shr 48);
  a_out[a_index + 7] := Byte(a_in shr 56);
end;

function ReverseBytesUInt32(value: UInt32): UInt32;
begin
result := (value and UInt32($000000FF)) shl 24 or (value and UInt32($0000FF00)
    ) shl 8 or (value and UInt32($00FF0000)) shr 8 or
    (value and UInt32($FF000000)) shr 24;
end;

procedure swap_copy_str_to_u32(src: Pointer; src_index: Int32;dest: Pointer; dest_index: Int32; length: Int32);
var
  lsrc, ldest, lend: PCardinal;
  lbsrc: PByte;
  lLength: Int32;
begin
  // if all pointers and length are 32-bits aligned
  if ((Int32(PByte(dest) - PByte(0)) or (PByte(src) - PByte(0)) or src_index or
    dest_index or length) and 3) = 0 then
  begin
    // copy memory as 32-bit words
    lsrc := PCardinal(PByte(src) + src_index);
    lend := PCardinal((PByte(src) + src_index) + length);
    ldest := PCardinal(PByte(dest) + dest_index);
    while lsrc <> lend do
    begin
      ldest^ := ReverseBytesUInt32(lsrc^);
      System.Inc(ldest);
      System.Inc(lsrc);
    end;
  end
  else
  begin
    lbsrc := (PByte(src) + src_index);
    lLength := length + dest_index;
    while dest_index < lLength do
    begin
      PByte(dest)[dest_index xor 3] := lbsrc^;
      System.Inc(lbsrc);
      System.Inc(dest_index);
    end;
  end;
end;


{ Digest256Bit }

function Digest256Bit.ToByteArray: TBytes;
begin
Setlength(Result,32);
swap_copy_str_to_u32(@Data[0],0,@Result[0],0,Length(Result));
end;

function Digest256Bit.ToString: AnsiString;
var t:Tbytes;
begin
Setlength(t,32);
swap_copy_str_to_u32(@Data[0],0,@t[0],0,Length(t));
Setlength(Result,length(t)*2);
BinToHex(t,Pansichar(Result),length(t));
end;


{ TSHA256 }

function TSHA256.Feed(a_data: PByte; a_length_a_data: Int32;var a_start_index, a_length: Int32; var a_processed_bytes: UInt64): Boolean;
var
  &Length: Int32;
begin
if (a_length_a_data = 0) then begin
                              result := false;
                              Exit;
                              end;
if (a_length = 0) then  begin
                        result := false;
                        Exit;
                        end;
Length := System.Length(Fm_buffer) - Fm_pos;
if (Length > a_length) then  Length := a_length;
System.Move(a_data[a_start_index], Fm_buffer[Fm_pos],Length * System.SizeOf(Byte));
Fm_pos := Fm_pos + Length;
a_start_index := a_start_index + Length;
a_length := a_length - Length;
a_processed_bytes := a_processed_bytes + UInt64(Length);
result := (Fm_pos = System.Length(Fm_buffer));
end;


function TSHA256.ComputeBytes(a_data: Tbytes): Digest256Bit;
begin
Initialize;
TransformBytes(a_data,0,length(a_data));
Finish;
Result:=Fm_State;
end;

constructor TSHA256.create;
begin
  Fm_block_size := 64;
  Fm_hash_size := 32;
  Fm_buffer_size := Int32(64 * 1024);
  System.SetLength(Fm_buffer, 64);
  System.FillChar(Fm_buffer[0], System.Length(Fm_buffer) * System.SizeOf(Byte), 0);
  Fm_pos := 0;
  Fm_processed_bytes := 0;
end;

procedure TSHA256.Initialize;
begin
  Fm_state.Data[0] := $6A09E667;
  Fm_state.Data[1] := $BB67AE85;
  Fm_state.Data[2] := $3C6EF372;
  Fm_state.Data[3] := $A54FF53A;
  Fm_state.Data[4] := $510E527F;
  Fm_state.Data[5] := $9B05688C;
  Fm_state.Data[6] := $1F83D9AB;
  Fm_state.Data[7] := $5BE0CD19;
  Fm_pos := 0;
  System.FillChar(Fm_buffer[0], System.Length(Fm_buffer) * System.SizeOf(Byte), 0);
  Fm_processed_bytes := 0;
end;

procedure TSHA256.TransformBytes(a_data: Tbytes;a_index, a_length: Int32);
var ptr_a_data: PByte;
begin
ptr_a_data := PByte(a_data);
if (Fm_Pos<>0) then begin
                    if (Feed(ptr_a_data, System.Length(a_data), a_index, a_length,Fm_processed_bytes)) then
                                  begin
                                  Fm_Pos:=0;
                                  TransformBlock(PByte(Fm_buffer), Length(Fm_buffer), 0);
                                  end;
                    end;
  while (a_length >= (Length(Fm_buffer))) do
    begin
    Fm_processed_bytes := Fm_processed_bytes + UInt64(Length(Fm_buffer));
    TransformBlock(ptr_a_data, Length(Fm_buffer), a_index);
    a_index := a_index + (Length(Fm_buffer));
    a_length := a_length - (LEngth(Fm_buffer));
    end;
if (a_length > 0) then Feed(ptr_a_data, System.Length(a_data), a_index, a_length,Fm_processed_bytes);
end;

procedure TSHA256.Finish;
var
  bits: UInt64;
  padindex: Int32;
  pad: Tbytes;
begin
bits := Fm_processed_bytes * 8;
if (Fm_pos < 56) then padindex := (56 - Fm_pos)
else padindex := (120 - Fm_pos);
System.SetLength(pad, padindex + 8);
pad[0] := $80;
bits := ReverseBytesUInt64(bits);
ReadUInt64AsBytesLE(bits, pad, padindex);
padindex := padindex + 8;
TransformBytes(pad, 0, padindex);
end;


procedure TSHA256.TransformBlock(a_data: PByte; a_data_length: Int32; a_index: Int32);
const
    s_K: array [0 .. 63] of UInt32 = ($428A2F98, $71374491, $B5C0FBCF,
      $E9B5DBA5, $3956C25B, $59F111F1, $923F82A4, $AB1C5ED5, $D807AA98,
      $12835B01, $243185BE, $550C7DC3, $72BE5D74, $80DEB1FE, $9BDC06A7,
      $C19BF174, $E49B69C1, $EFBE4786, $0FC19DC6, $240CA1CC, $2DE92C6F,
      $4A7484AA, $5CB0A9DC, $76F988DA, $983E5152, $A831C66D, $B00327C8,
      $BF597FC7, $C6E00BF3, $D5A79147, $06CA6351, $14292967, $27B70A85,
      $2E1B2138, $4D2C6DFC, $53380D13, $650A7354, $766A0ABB, $81C2C92E,
      $92722C85, $A2BFE8A1, $A81A664B, $C24B8B70, $C76C51A3, $D192E819,
      $D6990624, $F40E3585, $106AA070, $19A4C116, $1E376C08, $2748774C,
      $34B0BCB5, $391C0CB3, $4ED8AA4A, $5B9CCA4F, $682E6FF3, $748F82EE,
      $78A5636F, $84C87814, $8CC70208, $90BEFFFA, $A4506CEB, $BEF9A3F7,
      $C67178F2);
var
  A, B, C, D, E, F, G, H, T, T2: UInt32;
  r: Int32;
  data: array [0 .. 63] of UInt32;
begin
swap_copy_str_to_u32(a_data, a_index, @(data[0]), 0, 64);
A := Fm_state.Data[0];
B := Fm_state.Data[1];
C := Fm_state.Data[2];
D := Fm_state.Data[3];
E := Fm_state.Data[4];
F := Fm_state.Data[5];
G := Fm_state.Data[6];
H := Fm_state.Data[7];
   // Step 1
for r := 16 to 63 do   begin
                       T := data[r - 2];
                       T2 := data[r - 15];
                       data[r] := ((RotateRight32(T, 17)) xor (RotateRight32(T, 19))xor (T shr 10)) + data[r - 7] +
                       ((RotateRight32(T2, 7)) xor (RotateRight32(T2, 18))xor (T2 shr 3)) + data[r - 16];
                       end;
  // Step 2
for r := 0 to 63 do  begin
                     T := s_K[r] + data[r] + H +((RotateRight32(E, 6)) xor (RotateRight32(E, 11))
                            xor (RotateRight32(E, 25))) + ((E and F) xor (not E and G));
                            T2 := ((RotateRight32(A, 2)) xor (RotateRight32(A, 13))xor (RotateRight32(A, 22))) +
                            ((A and B) xor (A and C) xor (B and C));
                     H := G;
                     G := F;
                     F := E;
                     E := D + T;
                     D := C;
                     C := B;
                     B := A;
                     A := T + T2;
                     end;
Fm_state.Data[0] := Fm_state.Data[0] + A;
Fm_state.Data[1] := Fm_state.Data[1] + B;
Fm_state.Data[2] := Fm_state.Data[2] + C;
Fm_state.Data[3] := Fm_state.Data[3] + D;
Fm_state.Data[4] := Fm_state.Data[4] + E;
Fm_state.Data[5] := Fm_state.Data[5] + F;
Fm_state.Data[6] := Fm_state.Data[6] + G;
Fm_state.Data[7] := Fm_state.Data[7] + H;
System.FillChar(data, System.SizeOf(data), 0);
end;

function TSHA256.HashByteArray(B: Tbytes): Digest256Bit;
begin
ComputeBytes(B);
Move(Fm_state.Data[0],Result.Data[0],32);
end;

function TSHA256.HashFile(filename: String): Digest256Bit;
const
  BufSize = 64 * 1024;  // 64kb buffer
var
  fsFileToBeHashed: TFileStream;
  Buffer: Tbytes;
  IntFileSize : Uint64;
  i:integer;
begin
Setlength(Buffer,BufSize - 1);
try
    fsFileToBeHashed := TFileStream.Create(filename, fmOpenRead or fmShareDenyNone);
    IntFileSize      := fsFileToBeHashed.Size;
    Initialize;
  repeat
      i := fsFileToBeHashed.Read(Buffer, BufSize);
      if i <= 0 then break
          else TransformBytes(Buffer, 0,i);
  until false;
    Finish;
    Move(Fm_state.Data[0],Result.Data[0],32);
finally
    fsFileToBeHashed.free;
end;
end;

function TSHA256.HashString(S: String): Digest256Bit;
begin
ComputeBytes(TEncoding.UTF8.GetBytes(S));
Move(Fm_state.Data[0],Result.Data[0],32);
end;

function SHA256FileHash(_filename: string): TBytes;
var
   sha: TIdHashSHA256;
   fs: TFileStream;
begin
if TIdHashSHA256.IsAvailable then
    begin
     sha:= TIdHashSHA256.Create;
     try
      fs:= TFileStream.Create(_filename, fmOpenRead);
      try
       Result:= Tbytes(sha.HashStream(fs));
      finally
       sha.Free;
      end;
     finally
      fs.Free;
     end;
    end
else Result:=Sha256Hasher.HashFile(_filename).ToByteArray;
end;

function SHA256StringHash(_string: string): Tbytes;
var
   sha: TIdHashSHA256;
   z:Tbytes;
begin
if TIdHashSHA256.IsAvailable then
    begin
     sha:= TIdHashSHA256.Create;
     try
      z:= Tbytes(sha.HashString(_string));
      Result:=z;
     finally
      sha.Free;
     end;
    end
else Result:=Sha256Hasher.HashString(_string).ToByteArray;
end;

function SHA256BytesHash(_bytes: Tbytes): Tbytes;
var
   sha: TIdHashSHA256;
   var z:Tbytes;
begin
if TIdHashSHA256.IsAvailable then
    begin
     sha:= TIdHashSHA256.Create;
     try
      z:= Tbytes(sha.HashBytes (Tidbytes(_bytes)));
      result:=z;
     finally
      sha.Free;
     end;
    end
else Result:=Sha256Hasher.HashByteArray(_bytes).ToByteArray;
end;

begin
if not LoadOpenSSLLibrary then Sha256Hasher:=TSHA256.create;
end.
