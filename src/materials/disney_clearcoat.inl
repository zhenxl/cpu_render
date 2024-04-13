#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    // Homework 1: implement this!
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3d h = normalize(dir_in + dir_out);
    Real half_angle = abs(dot(h, dir_out));
    // Real out_angle  = abs(dot(frame.n, dir_out));
    Real in_angle   = abs(dot(frame.n, dir_in));
    Real R0 = Schlick_approximation(1.5);
    Real Fc = R0 + (1 - R0) * pow((1.0 - half_angle), 5.0);
    Real Dc = normal_distribution(clearcoat_gloss, h, frame);
    Real Gc = Mask_shadowing(dir_in, frame) * Mask_shadowing(dir_out, frame);
    Real f_clearcoat = Fc * Dc * Gc / 4 * in_angle;


    return make_const_spectrum(f_clearcoat);
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3d h = normalize(dir_in + dir_out);
    Real Dc = normal_distribution(clearcoat_gloss, h, frame);
    Real half_angle = abs(dot(h, dir_out));
    Real n_dot_h = abs(dot(h, frame.n));
    return Dc * n_dot_h / (4.0 * half_angle);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = (1.0 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
    Real alpha_2 = alpha * alpha;
    Real cos_h_elevation = sqrt((1 - pow(alpha_2, 1 - rnd_param_uv.x)) / (1 - alpha_2));
    Real h_azimuth = 2 * c_PI * rnd_param_uv.y;
    Real h_elevation = acos(cos_h_elevation);
    Vector3d h_l = Vector3d();
    h_l.x = sin(h_elevation) * cos(h_azimuth);
    h_l.y = sin(h_elevation) * sin(h_azimuth);
    h_l.z = cos_h_elevation;

    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, h_l);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, alpha /* roughness */
    };

    return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
