Add-Type @"
  using System;
  using System.Runtime.InteropServices;
  public class Win32Lib {
    [DllImport("kernel32.dll")]
    [return: MarshalAs(UnmanagedType.Bool)]
    public static extern bool IsProcessorFeaturePresent(uint ProcessorFeature);
  }
"@

enum ProcessorFeature {
    PF_SSE4_2_INSTRUCTIONS_AVAILABLE = 38
    PF_AVX_INSTRUCTIONS_AVAILABLE = 40
}

$supportedFeatures = @()

if([Win32Lib]::IsProcessorFeaturePresent(
    [uint32][ProcessorFeature]::PF_SSE4_2_INSTRUCTIONS_AVAILABLE
)) {
    $supportedFeatures += "sse4_2"
}

if([Win32Lib]::IsProcessorFeaturePresent(
    [uint32][ProcessorFeature]::PF_AVX_INSTRUCTIONS_AVAILABLE
)) {
    $supportedFeatures += "avx"
}

Write-Output ($supportedFeatures -join ";")
