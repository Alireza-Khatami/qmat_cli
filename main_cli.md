● You want: when you double-click a MAT vertex in the Polyscope window, the camera centers and zooms to focus on that
  point — the same as MeshLab's "Center on point" behavior when you double-click a surface point.

  In MeshLab double-click:
  1. Sets the trackball/pivot center to the clicked point (so rotation orbits around it)
  2. Smoothly flies the camera toward that point so it fills a reasonable portion of the view

  In Polyscope terms this maps to:
  - polyscope::view::lookAt(eye, target, up) — reposition the camera to look at the selected vertex
  - Or polyscope::view::resetCameraToHomeView() won't work since it's global

  Polyscope exposes:
  polyscope::view::lookAt(
      glm::vec3 cameraPos,   // where the camera sits
      glm::vec3 target,      // the point to look at (the MAT vertex)
      bool flyTo = false     // if true: animate; if false: snap
  );

  So the behavior would be:
  1. User double-clicks a MAT vertex (detected via Polyscope's pick system)
  2. Get that vertex's 3D position
  3. Call polyscope::view::lookAt(currentPos moved closer to target, target) — keeping the current viewing direction but
   recentering on the vertex and pulling the camera to a comfortable distance

  The tricky part: Polyscope doesn't natively expose "double-click" as a separate event from single-click. We'd need to
  detect it manually by timing two consecutive picks on the same vertex within ~300ms.




