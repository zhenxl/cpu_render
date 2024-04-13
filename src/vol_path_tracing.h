#pragma once

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if (vertex_) {
        PathVertex vertex = *vertex_;
        int current_medium_id = vertex.exterior_medium_id;
        auto medium = scene.media[current_medium_id];
        Spectrum sigma_a = get_sigma_a(medium, vertex.position);
        Real t = distance(ray.org, vertex.position);
        Spectrum transmittance = exp(-sigma_a * t);
        Spectrum emitted_light = make_zero_spectrum();
        if (is_light(scene.shapes[vertex.shape_id])) {
            emitted_light = emission(vertex, -ray.dir, scene);
        }
        return emitted_light * transmittance ;
    }
    return make_zero_spectrum();
}

Spectrum L_s1(const Ray& ray, const Medium& medium, const Light& light, const PointAndNormal& point_on_light, const Vector3& position, const Scene& scene, int light_id) {
    Vector3d dir_light = normalize(point_on_light.position - position);
    auto Geom = max(-dot(dir_light, point_on_light.normal), Real(0)) /
                        distance_squared(point_on_light.position, position);
    Ray shadow_ray{position, dir_light, 
                               get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(point_on_light.position, position)};
    Real visable = !occluded(scene, shadow_ray);
    Spectrum Le = emission(light, -dir_light, Real(0), point_on_light, scene);
    Spectrum sigma_t = get_sigma_a(medium, position) + get_sigma_s(medium, position);
    Spectrum exp_t = exp(-sigma_t * distance(point_on_light.position, position));
    PhaseFunction phase = get_phase_function(medium);
    Spectrum phase_part = eval(phase, -ray.dir, dir_light);
    Spectrum L_estimate = phase_part * Le * exp_t * Geom * visable;
    Real pdf = pdf_point_on_light(light, point_on_light, position, scene) * light_pmf(scene, light_id);
    return L_estimate / pdf;
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if (vertex_) {
        PathVertex vertex = *vertex_;
        int current_medium_id = vertex.exterior_medium_id;
        auto medium = scene.media[current_medium_id];
        Spectrum sigma_a = get_sigma_a(medium, vertex.position);
        Spectrum sigma_s = get_sigma_s(medium, vertex.position);
        Real sigma_t = (sigma_a + sigma_s).x;
        Real u = next_pcg32_real<Real>(rng);
        Real t = -log(1.0 - u) / sigma_t;
        Real t_hit = distance(ray.org, vertex.position);
        if (t < t_hit) {
            Real trans_pdf = exp(-sigma_t * t) * sigma_t;
            Real transmittance =  exp(-sigma_t * t);

            Vector3d p = ray.org + ray.dir * t;
            int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
            const Light &light = scene.lights[light_id];
            Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real shape_w = next_pcg32_real<Real>(rng);
            PointAndNormal point_on_light =
                sample_point_on_light(light, vertex.position, light_uv, shape_w, scene);
            auto L = L_s1(ray, medium, light, point_on_light, p, scene, light_id);
            return (transmittance / trans_pdf) * sigma_s * L;
        } else {
            Real trans_pdf = exp(-sigma_t * t_hit) ;
            Real transmittance =  exp(-sigma_t * t_hit);
            Spectrum Le = make_zero_spectrum();
            if (is_light(scene.shapes[vertex.shape_id])) {
                Le = emission(vertex, -ray.dir, scene);
            }
             return (transmittance / trans_pdf) * Le;
            //  return make_zero_spectrum();
        }
    } else {
        Real u = next_pcg32_real<Real>(rng);
        int current_medium_id = scene.camera.medium_id;
        auto medium = scene.media[current_medium_id];
        Spectrum sigma_a = get_sigma_a(medium, Vector3d{1,2,3});
        Spectrum sigma_s = get_sigma_s(medium, Vector3d{4,5,6});
        Real sigma_t = (sigma_a + sigma_s).x;
        Real t = -log(1.0 - u) / sigma_t;
        Real trans_pdf = exp(-sigma_t * t) * sigma_t;
        Real transmittance =  exp(-sigma_t * t);

        Vector3d p = ray.org + ray.dir * t;
        int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
        const Light &light = scene.lights[light_id];
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real shape_w = next_pcg32_real<Real>(rng);
        PointAndNormal point_on_light =
            sample_point_on_light(light, p, light_uv, shape_w, scene);
        auto L = L_s1(ray, medium, light, point_on_light, p, scene, light_id);
        return (transmittance / trans_pdf) * sigma_s * L;
    }
}

int update_medium(PathVertex vertex, Ray ray, int medium_id) {
    int medium = medium_id;
    if (vertex.interior_medium_id != vertex.exterior_medium_id) {
        if (dot(ray.dir, vertex.geometric_normal) > 0) {
            medium = vertex.exterior_medium_id;
        } else {
            medium = vertex.interior_medium_id;
        }
    }
    return medium;
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int current_medium_id = scene.camera.medium_id;
    // auto current_medium = scene.media[medium_id];
    Spectrum current_path_throughput = make_const_spectrum(1.0);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    int max_depth = scene.options.max_depth;
    while (true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Real trans_pdf = 1.0;
        Real transmittance = 1.0;
        if (current_medium_id != -1) {
            auto medium = scene.media[current_medium_id];
            Spectrum sigma_a = get_sigma_a(medium, Vector3d{1, 2, 3});
            Spectrum sigma_s = get_sigma_s(medium, Vector3d{1, 2, 3});
            Real sigma_t = (sigma_a + sigma_s).x;
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1.0 - u) / sigma_t;
            if (vertex_) {
                PathVertex vertex = *vertex_;
                Real t_hit = distance(ray.org, vertex.position);
                if (t < t_hit) {
                    scatter = true;
                    trans_pdf = exp(-sigma_t * t) * sigma_t;
                    transmittance = exp(-sigma_t * t);
                    ray.org = ray.org + t * ray.dir;
                } else {
                    scatter = false;
                    trans_pdf = exp(-sigma_t * t_hit);
                    transmittance =  exp(-sigma_t * t_hit);
                    ray.org = ray.org + t_hit * ray.dir;
                }
            } else {
                scatter = true;
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t) ;
                ray.org = ray.org + t * ray.dir;
            }
        }
        current_path_throughput *= (transmittance / trans_pdf);

        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (is_light(scene.shapes[vertex.shape_id])) {
                auto Le = emission(vertex, -ray.dir, scene);
                radiance += current_path_throughput * Le;
            }
        }

        if (bounces == max_depth - 1  && max_depth != -1) {
            break;
        }

        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (vertex.material_id == -1) {
                current_medium_id = update_medium(vertex, ray, current_medium_id);
                ray.org = vertex.position + ray.dir * get_intersection_epsilon(scene);
                bounces += 1;
                continue;
            }
        }

        if (scatter) {
            auto medium = scene.media[current_medium_id];
            PhaseFunction phase = get_phase_function(medium);
            std::optional<Vector3> next_dir_ = sample_phase_function(phase, -ray.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)));
            if (next_dir_) {
                Real sigma_s = get_sigma_s(medium, ray.org).x;
                Vector3 next_dir = *next_dir_;
                current_path_throughput *= (eval(phase, -ray.dir, next_dir) / pdf_sample_phase(phase, -ray.dir, next_dir)) * sigma_s;
                ray = Ray{ray.org + next_dir * get_intersection_epsilon(scene), next_dir, Real(0), infinity<Real>()};
            }
        } else {
            break;
        }

        Real rr_prob = 1.0;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(current_path_throughput.x, 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                // Terminate the path
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }
    return radiance;
}

Spectrum next_event_estimation(Vector3 pos, const Ray& ray, int medium_id, int bounces, const Scene &scene, pcg32_state& rng) {
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    const Light &light = scene.lights[light_id];
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real shape_w = next_pcg32_real<Real>(rng);
    PointAndNormal p_prime =
        sample_point_on_light(light, pos, light_uv, shape_w, scene);
    Real T_light = 1.0;
    int shadow_medium_id = medium_id;
    int shadow_bounce = 0;
    Real p_trans_dir = 1.0;
    int max_depth = scene.options.max_depth;
    Vector3 orig_pos = pos;
    while (true) {
        Ray shadow_ray{pos, normalize(p_prime.position - pos),
                       get_shadow_epsilon(scene),
                        (1 - get_shadow_epsilon(scene)) *
                        distance(p_prime.position, pos)};
        RayDifferential ray_diff{0.0, 0.0};
        std::optional<PathVertex> vertex_ = intersect(scene, shadow_ray, ray_diff);
        Real next_t = distance(pos, p_prime.position);
        if (vertex_) {
            PathVertex v = *vertex_;
            next_t = distance(pos, v.position);
        }

        if (shadow_medium_id != -1) {
            const Medium& medium = scene.media[shadow_medium_id];
            Real sigma_a = get_sigma_a(medium, Vector3{1, 2, 3}).x;
            Real sigma_s = get_sigma_s(medium, Vector3{1, 2, 3}).x;
            Real sigma_t = sigma_a + sigma_s;
            T_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t * next_t);
        }

        if (!vertex_) {
            break;
        } else {
            PathVertex vertex = *vertex_;
            if (vertex.material_id >= 0) {
                return make_zero_spectrum();
            }
        }

        shadow_bounce ++;

        if (max_depth != -1  && bounces + shadow_bounce + 1 >= max_depth) {
            return make_zero_spectrum();
        }

        shadow_medium_id = update_medium(*vertex_, shadow_ray, shadow_medium_id);
        pos = pos + next_t * shadow_ray.dir;
    }


    if (T_light > 0) {
        Vector3d dir_light = normalize(p_prime.position - orig_pos);
        auto Geom = max(-dot(dir_light, p_prime.normal), Real(0)) /
                        distance_squared(p_prime.position, pos);
        Medium medium = scene.media[medium_id];
        PhaseFunction phase = get_phase_function(medium);
        Spectrum f = eval(phase, -ray.dir, dir_light);
        Spectrum Le = emission(light, -dir_light, Real(0), p_prime, scene);
        
        Real pdf_nee = pdf_point_on_light(light, p_prime, orig_pos, scene) * light_pmf(scene, light_id);

        Spectrum contrib = T_light * Geom * f * Le / pdf_nee;
        Real pdf_phase = pdf_sample_phase(phase, -ray.dir, -dir_light) * Geom * p_trans_dir;
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
        return w * contrib ;
    }
    return make_zero_spectrum();
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int current_medium_id = scene.camera.medium_id;
    // auto current_medium = scene.media[medium_id];
    Spectrum current_path_throughput = make_const_spectrum(1.0);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    int max_depth = scene.options.max_depth;

    bool never_scatter = true;
    Real dir_pdf = 0.0;
    Vector3 nee_p_cache;
    Real multi_trans_pdf = 1.0;
    while (true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Real trans_pdf = 1.0;
        Real transmittance = 1.0;
        if (current_medium_id != -1) {
            auto medium = scene.media[current_medium_id];
            Spectrum sigma_a = get_sigma_a(medium, Vector3d{1, 2, 3});
            Spectrum sigma_s = get_sigma_s(medium, Vector3d{1, 2, 3});
            Real sigma_t = (sigma_a + sigma_s).x;
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1.0 - u) / sigma_t;
            if (vertex_) {
                PathVertex vertex = *vertex_;
                Real t_hit = distance(ray.org, vertex.position);
                if (t < t_hit) {
                    scatter = true;
                    trans_pdf = exp(-sigma_t * t) * sigma_t;
                    transmittance = exp(-sigma_t * t);
                    ray.org = ray.org + t * ray.dir;
                } else {
                    scatter = false;
                    trans_pdf = exp(-sigma_t * t_hit);
                    transmittance =  exp(-sigma_t * t_hit);
                    ray.org = ray.org + t_hit * ray.dir;
                } 
            } else {
                scatter = true;
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);
                ray.org = ray.org + t * ray.dir;
            }
        }
        current_path_throughput *= (transmittance / trans_pdf);

        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (never_scatter) {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    auto Le = emission(vertex, -ray.dir, scene);
                    radiance += current_path_throughput * Le;
                }
            } else {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    const Light &light = scene.lights[light_id];
                    PointAndNormal light_point{vertex.position, vertex.geometric_normal};
                    Real pdf_nee = pdf_point_on_light(light, light_point, nee_p_cache, scene) * light_pmf(scene, light_id);
                    Vector3 omega_prime = normalize(vertex.position - nee_p_cache);
                    Real top = abs(dot(omega_prime, vertex.geometric_normal));
                    Real bottom = length_squared(nee_p_cache - vertex.position);
                    Real G = top / bottom;
                    
                    Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    
                    Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * w;
                }
            }
        }

        if (bounces == max_depth - 1  && max_depth != -1) {
            break;
        }

        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (vertex.material_id == -1) {
                current_medium_id = update_medium(vertex, ray, current_medium_id);
                bounces += 1;
                ray.org = vertex.position + ray.dir * get_intersection_epsilon(scene);
                multi_trans_pdf *= trans_pdf;
                continue;
            }
        }

        if (scatter) {
            never_scatter = false;
            auto medium = scene.media[current_medium_id];
            PhaseFunction phase = get_phase_function(medium);
            std::optional<Vector3> next_dir_ = sample_phase_function(phase, -ray.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)));
            if (next_dir_) {
                Real sigma_s = get_sigma_s(medium, ray.org).x;
                Spectrum nee = next_event_estimation(ray.org, ray, current_medium_id, bounces, scene, rng);
                radiance += current_path_throughput * nee * sigma_s;
                Vector3 next_dir = *next_dir_;
                dir_pdf =  pdf_sample_phase(phase, ray.dir, next_dir);
                current_path_throughput *= ((eval(phase, -ray.dir, next_dir)) / dir_pdf) * sigma_s;
                nee_p_cache = ray.org;
                multi_trans_pdf = Real(1);
                ray = Ray{ray.org + next_dir * get_intersection_epsilon(scene), next_dir, Real(0), infinity<Real>()};
                
            }
        } else {
            //no scatter, we hit something, now need to brdf
            if (vertex_) {
                PathVertex vertex = *vertex_;

                const Material &mat = scene.materials[vertex.material_id];
                Vector3 dir_view = -ray.dir;
                Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
                std::optional<BSDFSampleRecord> bsdf_sample_ =
                        sample_bsdf(mat,
                                    dir_view,
                                    vertex,
                                    scene.texture_pool,
                                    bsdf_rnd_param_uv,
                                    bsdf_rnd_param_w);
                if (!bsdf_sample_) {
                        // BSDF sampling failed. Abort the loop.
                        break;
                }
                const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
                Vector3 dir_bsdf = bsdf_sample.dir_out;
                ray = Ray{vertex.position, dir_bsdf, get_intersection_epsilon(scene), infinity<Real>()};
                Spectrum f = eval(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
                Real pdf = pdf_sample_bsdf(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
                current_path_throughput *= f / pdf;
                dir_pdf = pdf;
                nee_p_cache = ray.org;
                multi_trans_pdf = Real(1);
            }
        }

        Real rr_prob = 1.0;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(current_path_throughput.x, 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                // Terminate the path
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }
    return radiance;
}

Spectrum next_event_estimation_brdf(Vector3 pos, const Ray& ray, int medium_id, int bounces, const Scene &scene, pcg32_state& rng, PathVertex& v) {
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    const Light &light = scene.lights[light_id];
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real shape_w = next_pcg32_real<Real>(rng);
    PointAndNormal p_prime =
        sample_point_on_light(light, pos, light_uv, shape_w, scene);
    Real T_light = 1.0;
    int shadow_medium_id = medium_id;
    int shadow_bounce = 0;
    Real p_trans_dir = 1.0;
    int max_depth = scene.options.max_depth;
    Vector3 orig_pos = pos;
    while (true) {
        Ray shadow_ray{pos, normalize(p_prime.position - pos),
                       get_shadow_epsilon(scene),
                        (1 - get_shadow_epsilon(scene)) *
                        distance(p_prime.position, pos)};
        RayDifferential ray_diff{0.0, 0.0};
        std::optional<PathVertex> vertex_ = intersect(scene, shadow_ray, ray_diff);
        Real next_t = distance(pos, p_prime.position);
        if (vertex_) {
            PathVertex v = *vertex_;
            next_t = distance(pos, v.position);
        }

        if (shadow_medium_id != -1) {
            const Medium& medium = scene.media[shadow_medium_id];
            Real sigma_a = get_sigma_a(medium, Vector3{1, 2, 3}).x;
            Real sigma_s = get_sigma_s(medium, Vector3{1, 2, 3}).x;
            Real sigma_t = sigma_a + sigma_s;
            T_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t * next_t);
        }

        if (!vertex_) {
            break;
        } else {
            PathVertex vertex = *vertex_;
            if (vertex.material_id >= 0) {
                return make_zero_spectrum();
            }
        }

        shadow_bounce ++;

        if (max_depth != -1  && bounces + shadow_bounce + 1 >= max_depth) {
            return make_zero_spectrum();
        }

        shadow_medium_id = update_medium(*vertex_, shadow_ray, shadow_medium_id);
        pos = pos + next_t * shadow_ray.dir;
    }


    if (T_light > 0) {
        Vector3d dir_light = normalize(p_prime.position - orig_pos);
        auto Geom = max(-dot(dir_light, p_prime.normal), Real(0)) /
                        distance_squared(p_prime.position, pos);
        Medium medium = scene.media[medium_id];
        PhaseFunction phase = get_phase_function(medium);
        const Material& mat = scene.materials[v.material_id];
        Spectrum f = eval(mat,  -ray.dir, dir_light, v, scene.texture_pool);
        Spectrum Le = emission(light, -dir_light, Real(0), p_prime, scene);
        
        Real pdf_nee = pdf_point_on_light(light, p_prime, orig_pos, scene) * light_pmf(scene, light_id);

        Spectrum contrib = T_light * Geom * f * Le / pdf_nee;
        // Real pdf_phase = pdf_sample_phase(phase, -ray.dir, -dir_light) * Geom * p_trans_dir;
       
        Real pdf_brdf = pdf_sample_bsdf(mat, -ray.dir, dir_light, v, scene.texture_pool) * Geom * p_trans_dir;;
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_brdf * pdf_brdf);
        return w * contrib;
    }
    return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
     // Homework 2: implememt this!
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int current_medium_id = scene.camera.medium_id;
    // auto current_medium = scene.media[medium_id];
    Spectrum current_path_throughput = make_const_spectrum(1.0);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    int max_depth = scene.options.max_depth;

    bool never_scatter = true;
    Real dir_pdf = 0.0;
    Vector3 nee_p_cache;
    Real multi_trans_pdf = 1.0;
    while (true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Real trans_pdf = 1.0;
        Real transmittance = 1.0;
        if (current_medium_id != -1) {
            auto medium = scene.media[current_medium_id];
            Spectrum sigma_a = get_sigma_a(medium, Vector3d{1, 2, 3});
            Spectrum sigma_s = get_sigma_s(medium, Vector3d{1, 2, 3});
            Real sigma_t = (sigma_a + sigma_s).x;
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1.0 - u) / sigma_t;
            if (vertex_) {
                PathVertex vertex = *vertex_;
                Real t_hit = distance(ray.org, vertex.position);
                if (t < t_hit) {
                    scatter = true;
                    trans_pdf = exp(-sigma_t * t) * sigma_t;
                    transmittance = exp(-sigma_t * t);
                    ray.org = ray.org + t * ray.dir;
                } else {
                    scatter = false;
                    trans_pdf = exp(-sigma_t * t_hit);
                    transmittance =  exp(-sigma_t * t_hit);
                    ray.org = ray.org + t_hit * ray.dir;
                } 
            } else {
                scatter = true;
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);
                ray.org = ray.org + t * ray.dir;
            }
        }
        current_path_throughput *= (transmittance / trans_pdf);

        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (never_scatter) {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    auto Le = emission(vertex, -ray.dir, scene);
                    radiance += current_path_throughput * Le;
                }
            } else {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    const Light &light = scene.lights[light_id];
                    PointAndNormal light_point{vertex.position, vertex.geometric_normal};
                    Real pdf_nee = pdf_point_on_light(light, light_point, nee_p_cache, scene) * light_pmf(scene, light_id);
                    Vector3 omega_prime = normalize(vertex.position - nee_p_cache);
                    Real top = abs(dot(omega_prime, vertex.geometric_normal));
                    Real bottom = length_squared(nee_p_cache - vertex.position);
                    Real G = top / bottom;
                    
                    Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    
                    Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * w;
                }
            }
        }

        if (bounces == max_depth - 1  && max_depth != -1) {
            break;
        }

        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (vertex.material_id == -1) {
                current_medium_id = update_medium(vertex, ray, current_medium_id);
                bounces += 1;
                ray.org = vertex.position + ray.dir * get_intersection_epsilon(scene);
                multi_trans_pdf *= trans_pdf;
                continue;
            }
        }

        if (scatter) {
            never_scatter = false;
            auto medium = scene.media[current_medium_id];
            PhaseFunction phase = get_phase_function(medium);
            std::optional<Vector3> next_dir_ = sample_phase_function(phase, -ray.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)));
            if (next_dir_) {
                Real sigma_s = get_sigma_s(medium, ray.org).x;
                Spectrum nee = next_event_estimation(ray.org, ray, current_medium_id, bounces, scene, rng);
                radiance += current_path_throughput * nee * sigma_s;
                Vector3 next_dir = *next_dir_;
                dir_pdf =  pdf_sample_phase(phase, ray.dir, next_dir);
                current_path_throughput *= ((eval(phase, -ray.dir, next_dir)) / dir_pdf) * sigma_s;
                nee_p_cache = ray.org;
                multi_trans_pdf = Real(1);
                ray = Ray{ray.org + next_dir * get_intersection_epsilon(scene), next_dir, Real(0), infinity<Real>()};
                
            }
        } else {
            //no scatter, we hit something, now need to brdf
            if (vertex_) {
                PathVertex vertex = *vertex_;

                //nee 
                Spectrum nee_brdf = next_event_estimation_brdf(ray.org, ray, current_medium_id, bounces, scene, rng, vertex);
                radiance += current_path_throughput * nee_brdf;
                never_scatter = false;
                const Material &mat = scene.materials[vertex.material_id];
                Vector3 dir_view = -ray.dir;
                Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
                std::optional<BSDFSampleRecord> bsdf_sample_ =
                        sample_bsdf(mat,
                                    dir_view,
                                    vertex,
                                    scene.texture_pool,
                                    bsdf_rnd_param_uv,
                                    bsdf_rnd_param_w);
                if (!bsdf_sample_) {
                        // BSDF sampling failed. Abort the loop.
                        break;
                }
                const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
                Vector3 dir_bsdf = bsdf_sample.dir_out;
                ray = Ray{vertex.position, dir_bsdf, get_intersection_epsilon(scene), infinity<Real>()};
                Spectrum f = eval(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
                Real pdf = pdf_sample_bsdf(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
                current_path_throughput *= f / pdf;
                dir_pdf = pdf;
                nee_p_cache = ray.org;
                multi_trans_pdf = Real(1);
            }
        }

        Real rr_prob = 1.0;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(current_path_throughput.x, 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                // Terminate the path
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }
    return radiance;
}

Spectrum next_event_estimation_final(Vector3 pos, const Ray& ray, int medium_id, int bounces, const Scene &scene, pcg32_state& rng) {
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    const Light &light = scene.lights[light_id];
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real shape_w = next_pcg32_real<Real>(rng);
    PointAndNormal p_prime =
        sample_point_on_light(light, pos, light_uv, shape_w, scene);
    Spectrum T_light = make_const_spectrum(1.0);
    int shadow_medium_id = medium_id;
    int shadow_bounce = 0;
    Spectrum p_trans_dir = make_const_spectrum(1.0);
    Spectrum p_trans_nee = make_const_spectrum(1.0);
    int max_depth = scene.options.max_depth;
    Vector3 orig_pos = pos;
    while (true) {
        Ray shadow_ray{pos, normalize(p_prime.position - pos),
                       get_shadow_epsilon(scene),
                        (1 - get_shadow_epsilon(scene)) *
                        distance(p_prime.position, pos)};
        RayDifferential ray_diff{0.0, 0.0};
        std::optional<PathVertex> vertex_ = intersect(scene, shadow_ray, ray_diff);
        Real next_t = distance(pos, p_prime.position);
        if (vertex_) {
            PathVertex v = *vertex_;
            next_t = distance(pos, v.position);
        }

        if (shadow_medium_id != -1) {
            int channel = std::clamp(int(next_pcg32_real<Real>(rng))*3, 0, 2);
            int iteration = 0;
            Real accum_t = 0.0;
            const Medium& medium = scene.media[shadow_medium_id];
            Spectrum majorant = get_majorant(medium, shadow_ray);

            while (true) {
                if (majorant[channel] <= 0) {
                    break;
                } 

                if (iteration >= scene.options.max_null_collisions) {
                    break;
                }

                Real t = -log(1.0 - next_pcg32_real<Real>(rng)) / majorant[channel];

                Real dt = next_t - accum_t;
                accum_t = std::min(accum_t + t, next_t);
                
                Vector3 curr_point = ray.org + accum_t * ray.dir;
                Spectrum sigma_a = get_sigma_a(medium, curr_point);
                Spectrum sigma_s = get_sigma_s(medium, curr_point);
                Spectrum sigma_t = sigma_a + sigma_s; 
            }

            Real sigma_a = get_sigma_a(medium, Vector3{1, 2, 3}).x;
            Real sigma_s = get_sigma_s(medium, Vector3{1, 2, 3}).x;
            Real sigma_t = sigma_a + sigma_s;
            T_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t * next_t);
        }

        if (!vertex_) {
            break;
        } else {
            PathVertex vertex = *vertex_;
            if (vertex.material_id >= 0) {
                return make_zero_spectrum();
            }
        }

        shadow_bounce ++;

        if (max_depth != -1  && bounces + shadow_bounce + 1 >= max_depth) {
            return make_zero_spectrum();
        }

        shadow_medium_id = update_medium(*vertex_, shadow_ray, shadow_medium_id);
        pos = pos + next_t * shadow_ray.dir;
    }


    /*if (T_light > 0) {
        Vector3d dir_light = normalize(p_prime.position - orig_pos);
        auto Geom = max(-dot(dir_light, p_prime.normal), Real(0)) /
                        distance_squared(p_prime.position, pos);
        Medium medium = scene.media[medium_id];
        PhaseFunction phase = get_phase_function(medium);
        Spectrum f = eval(phase, -ray.dir, dir_light);
        Spectrum Le = emission(light, -dir_light, Real(0), p_prime, scene);
        
        Real pdf_nee = pdf_point_on_light(light, p_prime, orig_pos, scene) * light_pmf(scene, light_id);

        Spectrum contrib = T_light * Geom * f * Le / pdf_nee;
        Real pdf_phase = pdf_sample_phase(phase, -ray.dir, -dir_light) * Geom * p_trans_dir;
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
        return w * contrib ;
    }*/
    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int current_medium_id = scene.camera.medium_id;
    // auto current_medium = scene.media[medium_id];
    Spectrum current_path_throughput = make_const_spectrum(1.0);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    int max_depth = scene.options.max_depth;

    bool never_scatter = true;
    Spectrum dir_pdf = make_zero_spectrum();
    Vector3 nee_p_cache;
    Spectrum multi_trans_pdf = make_const_spectrum(1.0);
    while (true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Spectrum trans_dir_pdf = make_const_spectrum(1.0);
        Spectrum transmittance = make_const_spectrum(1.0);
        Spectrum trans_nee_pdf = make_const_spectrum(1.0);
        Real t = 0.0;
        if (current_medium_id != -1) {
            auto medium = scene.media[current_medium_id];
            Spectrum majorant = get_majorant(medium, ray);
            int channel = std::clamp(int(next_pcg32_real<Real>(rng))*3, 0, 2);
            Real accum_t = 0.0;
            int iteration = 0;

            Real t_hit = infinity<Real>();

            if (vertex_) {
                PathVertex vertex = *vertex_;
                t_hit = distance(ray.org, vertex.position);
            }

            while (true) {
                if (majorant[channel] <= 0) {
                    break;
                }

                if (iteration >= scene.options.max_null_collisions) {
                    break;
                }

                t = -log(1.0 - next_pcg32_real<Real>(rng)) / majorant[channel];

                Real dt = t_hit - accum_t;
                accum_t = std::min(accum_t + t, t_hit);
                
                Vector3 curr_point = ray.org + accum_t * ray.dir;
                Spectrum sigma_a = get_sigma_a(medium, curr_point);
                Spectrum sigma_s = get_sigma_s(medium, curr_point);
                Spectrum sigma_t = sigma_a + sigma_s; 

                if (t < dt) {
                    Spectrum real_prob = sigma_t / majorant;

                    if (next_pcg32_real<Real>(rng) < real_prob[channel]) {
                        scatter = true;
                        never_scatter = false;
                        transmittance *= exp(-majorant * t)/ max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);
                        break;
                    } else {
                        Spectrum sigma_n = majorant - sigma_t;
                        transmittance *= exp(-majorant * t) * sigma_n / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);
                        trans_nee_pdf *= exp(-majorant * t) * majorant / max(majorant);
                    }
                } else {
                    // reach the surface
                    transmittance *= exp(-majorant * dt);
                    trans_dir_pdf *= exp(-majorant * dt);
                    trans_nee_pdf *= exp(-majorant * dt);
                    break;
                }
                iteration++;
            }
             ray.org += accum_t * ray.dir;
        }
        current_path_throughput *= (transmittance / trans_dir_pdf);

        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (never_scatter) {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    auto Le = emission(vertex, -ray.dir, scene);
                    radiance += current_path_throughput * Le;
                }
            } else {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    const Light &light = scene.lights[light_id];
                    PointAndNormal light_point{vertex.position, vertex.geometric_normal};
                    Spectrum pdf_nee = pdf_point_on_light(light, light_point, nee_p_cache, scene) * light_pmf(scene, light_id) * trans_nee_pdf;
                    Vector3 omega_prime = normalize(vertex.position - nee_p_cache);
                    Real top = abs(dot(omega_prime, vertex.geometric_normal));
                    Real bottom = length_squared(nee_p_cache - vertex.position);
                    Real G = top / bottom;
                    
                    Spectrum dir_pdf_ = trans_dir_pdf * multi_trans_pdf * G;
                    
                    Spectrum w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * w;
                }
            }
        }

        if (bounces == max_depth - 1  && max_depth != -1) {
            break;
        }

        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (vertex.material_id == -1) {
                current_medium_id = update_medium(vertex, ray, current_medium_id);
                bounces += 1;
                ray.org = vertex.position + ray.dir * get_intersection_epsilon(scene);
                multi_trans_pdf *= trans_dir_pdf;
                continue;
            }
        }

        if (scatter) {
            never_scatter = false;
            auto medium = scene.media[current_medium_id];
            PhaseFunction phase = get_phase_function(medium);
            std::optional<Vector3> next_dir_ = sample_phase_function(phase, -ray.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)));
            if (next_dir_) {
                Real sigma_s = get_sigma_s(medium, ray.org).x;
                Spectrum nee = next_event_estimation(ray.org, ray, current_medium_id, bounces, scene, rng);
                radiance += current_path_throughput * nee * sigma_s;
                Vector3 next_dir = *next_dir_;
                trans_dir_pdf =  make_const_spectrum(pdf_sample_phase(phase, ray.dir, next_dir));
                current_path_throughput *= ((eval(phase, -ray.dir, next_dir)) / trans_dir_pdf) * sigma_s;
                nee_p_cache = ray.org;
                multi_trans_pdf = make_const_spectrum(1.0);
                ray = Ray{ray.org + next_dir * get_intersection_epsilon(scene), next_dir, Real(0), infinity<Real>()};
                
            }
        } else {
            //no scatter, we hit something, now need to brdf
            if (vertex_) {
                PathVertex vertex = *vertex_;

                //nee 
                Spectrum nee_brdf = next_event_estimation_brdf(ray.org, ray, current_medium_id, bounces, scene, rng, vertex);
                radiance += current_path_throughput * nee_brdf;
                never_scatter = false;
                const Material &mat = scene.materials[vertex.material_id];
                Vector3 dir_view = -ray.dir;
                Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
                std::optional<BSDFSampleRecord> bsdf_sample_ =
                        sample_bsdf(mat,
                                    dir_view,
                                    vertex,
                                    scene.texture_pool,
                                    bsdf_rnd_param_uv,
                                    bsdf_rnd_param_w);
                if (!bsdf_sample_) {
                        // BSDF sampling failed. Abort the loop.
                        break;
                }
                const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
                Vector3 dir_bsdf = bsdf_sample.dir_out;
                ray = Ray{vertex.position, dir_bsdf, get_intersection_epsilon(scene), infinity<Real>()};
                Spectrum f = eval(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
                Real pdf = pdf_sample_bsdf(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
                current_path_throughput *= f / pdf;
            /*    dir_pdf = pdf;
                nee_p_cache = ray.org;
                multi_trans_pdf = Real(1);*/
            }
        }

        Real rr_prob = 1.0;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(current_path_throughput.x, 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                // Terminate the path
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }
    return radiance;
}
