
#ifndef _COMPRESSIONDECOMPRESSION_H_
#define _COMPRESSIONDECOMPRESSION_H_

#include <zfp.h>
#include <vector>

zfp_field* compressMatrix(zfp_stream* zfp, zfp_field* field, size_t bufsize, std::vector<double>& data);
zfp_field* decompressMatrix(zfp_stream* zfp, zfp_field* field, size_t bufsize, zfp_field* dec_field);

#endif // _COMPRESSIONDECOMPRESSION_H_
