'correct the typedef MC_HANDLE for 32-bit architecture

set objWS = CreateObject("Wscript.Shell")
Set fsob=CreateObject("Scripting.FileSystemObject")

old_text="typedef unsigned long long int MC_HANDLE"
new_text="typedef unsigned int  MC_HANDLE;"

strFileName = "../bin/mconf_matlab32.h"

Const FOR_READING = 1
Const FOR_WRITING = 2
strCheckForString = UCase("typedef unsigned long long int MC_HANDLE")
Set objFS = CreateObject("Scripting.FileSystemObject")
Set objTS = objFS.OpenTextFile(strFileName, FOR_READING)
strContents = objTS.ReadAll
objTS.Close
arrLines = Split(strContents, vbNewLine)
Set objTS = objFS.OpenTextFile(strFileName, FOR_WRITING)
For Each strLine In arrLines
  If (Left(UCase(LTrim(strLine)),Len(strCheckForString)) = strCheckForString) Then
    objTS.WriteLine new_text
  else
    objTS.WriteLine strLine
  End If
Next
