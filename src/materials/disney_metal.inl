#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
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
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 h = normalize(dir_in + dir_out);
    Real in_angle   = abs(dot(frame.n, dir_in));
    Real half_angle = abs(dot(h, dir_out));
    auto F_m  = base_color  + (1 - base_color) * pow(1 - half_angle, 5.0);
   
    Real Dm = Disney_GGX(anisotropic, h, roughness, frame);

    Real G_win  = Disney_Smith(anisotropic, dir_in, roughness, frame);
    Real G_wout = Disney_Smith(anisotropic, dir_out, roughness, frame);
    Real Gm = G_win * G_wout;

    auto f_metal = F_m * Gm * Dm / (4 * in_angle);

    return f_metal;
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 h = normalize(dir_in + dir_out);
    Real in_angle   = abs(dot(frame.n, dir_in));
    Real Dm = Disney_GGX(anisotropic, h, roughness, frame);

    Real G_win  = Disney_Smith(anisotropic, dir_in, roughness, frame);

    return Dm * G_win / (4 * in_angle);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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
      // Convert the incoming direction to local coordinates
        Vector3 local_dir_in = to_local(frame, dir_in);
        Real roughness = eval(
            bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
        // Clamp roughness to avoid numerical issues.
        roughness = std::clamp(roughness, Real(0.01), Real(1));
        // Real alpha = roughness * roughness;
        Real aspect = sqrt(1.0 - 0.9 * anisotropic);
        Real alpha_min = 0.0001;
        Real alpha_x = std::max(alpha_min, roughness*roughness / aspect);
        Real alpha_y = std::max(alpha_min, roughness * roughness * aspect);
        Vector3 local_micro_normal =
            sample_visible_normals(local_dir_in, alpha_x, alpha_y, rnd_param_uv);
        
        // Transform the micro normal to world space
        Vector3 half_vector = to_world(frame, local_micro_normal);
        // Reflect over the world space normal
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        return BSDFSampleRecord{
            reflected,
            Real(0) /* eta */, roughness /* roughness */
        };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
