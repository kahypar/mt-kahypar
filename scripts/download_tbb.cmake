include(FetchContent)
FetchContent_Populate(
  tbb
  URL https://github.com/oneapi-src/oneTBB/releases/download/v2020.2/tbb-2020.2-lin.tgz
  URL_HASH MD5=ffd526bcddbfbe342b841fe51a23ac86
  SOURCE_DIR external_tools/tbb
)