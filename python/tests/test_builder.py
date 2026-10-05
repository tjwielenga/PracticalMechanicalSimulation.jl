from __future__ import annotations

import tempfile
import tomllib
import unittest
from pathlib import Path

from practical_mechanical_simulation import Sim2D, Sim3D


class BuilderTests(unittest.TestCase):
    def test_planar_hierarchy_and_references(self):
        model = Sim2D.Model("test", title="Planar test")
        ground = Sim2D.ground(model, "ground")
        pin = Sim2D.marker(ground, "pin", position=[0, 0.0])
        body = Sim2D.rigid_body(
            model, "assembly.link", mass=1, inertia=0.2, position=[0.5, 0]
        )
        end = Sim2D.marker(body, "end", position=[-0.5, 0])
        joint = Sim2D.revolute(
            model,
            "assembly.joint",
            markers=[end, pin],
            rotation_coordinates=True,
        )
        Sim2D.state_selection(
            model, preferred_velocities=[Sim2D.variable(joint, "omega")]
        )
        doc = Sim2D.document(model)
        self.assertEqual(doc["assembly"]["link"]["end"]["type"], "marker")
        self.assertEqual(
            doc["assembly"]["joint"]["markers"],
            ["assembly.link.end", "ground.pin"],
        )
        self.assertEqual(
            doc["state_selection"]["preferred_velocities"],
            ["assembly.joint.omega"],
        )
        self.assertEqual(doc["assembly"]["link"]["position"], [0.5, 0.0])

    def test_spatial_named_builder_and_samples(self):
        model = Sim3D.Model("test")
        body = Sim3D.rigid_body(
            model,
            "body",
            mass=2.0,
            inertia=[1.0, 2.0, 3.0],
            position=[0.0, 0.0, 0.0],
        )
        Sim3D.graphics(body, shape="box", size=[1.0, 0.2, 0.1])
        Sim3D.simulation(
            model, start_time=1.0, end_time=3.0, frames_per_second=60
        )
        doc = Sim3D.document(model)
        self.assertEqual(doc["body"]["type"], "rigid_body")
        self.assertEqual(doc["body"]["graphics"]["shape"], "box")
        self.assertEqual(doc["simulation"]["output_samples"], 121)

    def test_toml_round_trip(self):
        model = Sim3D.Model("round_trip", title="Round trip")
        ground = Sim3D.ground(model, "ground")
        Sim3D.marker(ground, "origin", position=[0.0, 0.0, 0.0])
        Sim3D.analysis(model, mode="automatic")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.toml"
            Sim3D.write_model(path, model)
            parsed = tomllib.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(parsed, Sim3D.document(model))

    def test_cross_model_reference_is_rejected(self):
        first = Sim2D.Model("first")
        second = Sim2D.Model("second")
        first_ground = Sim2D.ground(first, "ground")
        second_ground = Sim2D.ground(second, "ground")
        foreign_marker = Sim2D.marker(first_ground, "pin")
        local_marker = Sim2D.marker(second_ground, "pin")
        with self.assertRaisesRegex(ValueError, "different model"):
            Sim2D.revolute(
                second, "joint", markers=[foreign_marker, local_marker]
            )

    def test_duplicate_and_invalid_marker_owner_are_rejected(self):
        model = Sim2D.Model("test")
        ground = Sim2D.ground(model, "ground")
        with self.assertRaisesRegex(ValueError, "already defined"):
            Sim2D.ground(model, "ground")
        joint = Sim2D.element(model, "revolute", "joint", markers=[])
        with self.assertRaisesRegex(ValueError, "must be owned"):
            Sim2D.marker(joint, "bad")
        self.assertEqual(str(ground), "ground")

    def test_beam_marker(self):
        model = Sim3D.Model("test")
        beam = Sim3D.flexible_beam(model, "beam", length=1.0)
        self.assertEqual(str(Sim3D.beam_marker(beam, "end_j")), "beam.end_j")
        self.assertEqual(str(Sim3D.marker(beam, "end_i")), "beam.end_i")
        self.assertNotIn("end_i", Sim3D.document(model)["beam"])
        with self.assertRaisesRegex(ValueError, "end_i"):
            Sim3D.beam_marker(beam, "quarter")
        with self.assertRaisesRegex(ValueError, "do not accept"):
            Sim3D.marker(beam, "cm", position=[0.0, 0.0, 0.0])


if __name__ == "__main__":
    unittest.main()
