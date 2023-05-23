// main.cpp
#include <iostream>
#include <zfp.h>

int main() {
    zfp_type type = zfp_type_double; // Set the data type.
    zfp_field* field = zfp_field_alloc(); // Create a ZFP field object.

    // Your ZFP library usage code here.

    zfp_field_free(field); // Don't forget to free the field object.

    std::cout << "ZFP library example" << std::endl;
    return 0;
}

