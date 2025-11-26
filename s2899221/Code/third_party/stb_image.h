// stb_image - public domain header
// https://github.com/nothings/stb
#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#endif
#ifndef STB_IMAGE_STATIC
#define STB_IMAGE_STATIC
#endif
#ifndef STBI_NO_STDIO
#define STBI_NO_STDIO
#endif
#define STBI_ONLY_JPEG
#define STBI_ONLY_PNG

#include <stddef.h>

#ifndef STBI_INCLUDE_STB_IMAGE_H
#define STBI_INCLUDE_STB_IMAGE_H
extern "C" {
unsigned char *stbi_load_from_memory(const unsigned char *buffer, int len, int *x, int *y, int *channels_in_file, int desired_channels);
void stbi_image_free(void *retval_from_stbi_load);
}
#endif


