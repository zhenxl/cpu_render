#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    Spectrum base_color = eval(bsdf.base_color, vertex.uv,  vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    auto lum_base_color = luminance(base_color);
    Spectrum C_tint;
    if (lum_base_color > 0) {
        C_tint = base_color / lum_base_color;
    } else {
        C_tint = make_const_spectrum(1.0);
    }
    Spectrum C_sheen = 1.0 - sheen_tint + sheen_tint * C_tint;
    Vector3 h = normalize(dir_in + dir_out);
    Real half_angle = abs(dot(h, dir_out));
    Real out_angle = abs(dot(frame.n, dir_out));
    Spectrum f_sheen = C_sheen * pow((1.0 - half_angle), 5.0) * out_angle;
    // Homework 1: implement this!
    return f_sheen;
}

Real pdf_sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(0) /* roughness */ };
}

TextureSpectrum get_texture_op::operator()(const DisneySheen &bsdf) const {
    return bsdf.base_color;
}
