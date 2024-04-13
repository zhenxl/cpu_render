Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 h = normalize(dir_in + dir_out);
    Real half_angle = abs(dot(h, dir_out));
    Real out_angle  = abs(dot(frame.n, dir_out));
    Real in_angle   = abs(dot(frame.n, dir_in));
    Real FD_90 = 0.5 + 2 * roughness * pow(half_angle, Real(2));
    Real FD_w_out = 1.0 + (FD_90 - 1.0) * pow((1.0 - out_angle), Real(5));
    Real FD_w_in  = 1.0 + (FD_90 - 1.0) * pow((1.0 - in_angle), Real(5));
    Spectrum f_base_diffuse = base_color  * c_INVPI * FD_w_in * FD_w_out * out_angle;

    Real FSS_90 = roughness * pow(out_angle, 2.0);
    Real FSS_w_out = 1.0 + (FSS_90 - 1.0) * pow((1.0 - out_angle), Real(5));
    Real FSS_w_in = 1.0 + (FSS_90 - 1.0) * pow((1.0 - in_angle), Real(5));
    Spectrum f_surface = 1.25 * base_color / c_PI * ((1.0 / (in_angle + out_angle) - 0.5) * FSS_w_in * FSS_w_out + 0.5) * out_angle;
    Real eval_subsurface = eval(bsdf.subsurface,  vertex.uv, vertex.uv_screen_size, texture_pool);

    Spectrum f_diffuse = eval_subsurface * f_surface + (1 - eval_subsurface) * f_base_diffuse;

    // Homework 1: implement this!
    return f_diffuse;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    Real eval_roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, eval_roughness  /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
