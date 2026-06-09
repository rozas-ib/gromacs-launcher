import os


def get_config_dir(cfg):
    config_dir = cfg.get("__config_dir")
    if not config_dir:
        raise ValueError("Configuration metadata is missing __config_dir")
    return config_dir


def resolve_project_path(cfg, path_value):
    if not isinstance(path_value, str) or not path_value.strip():
        raise ValueError("Expected a non-empty path string in the configuration")
    path_value = path_value.strip()
    if os.path.isabs(path_value):
        return path_value
    return os.path.abspath(os.path.join(get_config_dir(cfg), path_value))


def get_inputs_dir(cfg):
    return resolve_project_path(cfg, "inputs")


def get_scripts_dir(cfg):
    project_cfg = cfg.get("project_settings", {})
    return resolve_project_path(cfg, project_cfg.get("scripts_dir", "scripts"))


def get_template_dir(cfg):
    project_cfg = cfg.get("project_settings", {})
    return resolve_project_path(cfg, project_cfg.get("template_dir", "mdp_templates"))


def get_output_root(cfg):
    project_cfg = cfg.get("project_settings", {})
    output_root = project_cfg.get("output_root", "data")
    if not isinstance(output_root, str) or not output_root.strip():
        raise ValueError("project_settings.output_root must be a non-empty string")
    return resolve_project_path(cfg, output_root)


def get_case_dir(output_root, label):
    return os.path.join(output_root, label)


def get_replica_dir(case_dir, replica_idx):
    return os.path.join(case_dir, f"rep_{replica_idx}")


def get_stage_dir(replica_dir, stage_name):
    return os.path.join(replica_dir, stage_name)
