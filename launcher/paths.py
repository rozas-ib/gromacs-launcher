import os


def get_output_root(cfg):
    project_cfg = cfg.get("project_settings", {})
    output_root = project_cfg.get("output_root", "data")
    if not isinstance(output_root, str) or not output_root.strip():
        raise ValueError("project_settings.output_root must be a non-empty string")
    return os.path.abspath(output_root.strip())


def get_case_dir(output_root, label):
    return os.path.join(output_root, label)


def get_replica_dir(case_dir, replica_idx):
    return os.path.join(case_dir, f"rep_{replica_idx}")


def get_stage_dir(replica_dir, stage_name):
    return os.path.join(replica_dir, stage_name)

