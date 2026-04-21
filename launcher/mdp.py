def modify_mdp(template_path, output_path, replacements):
    """Replace MDP keys with tailored values."""
    with open(template_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    new_lines = []
    for line in lines:
        cleaned = line.split(";")[0].strip()
        if "=" in cleaned:
            key = cleaned.split("=")[0].strip()
            if key in replacements:
                new_lines.append(f"{key:<24}= {replacements[key]} ; Tailored\n")
                continue
        new_lines.append(line)
    with open(output_path, "w", encoding="utf-8") as f:
        f.writelines(new_lines)

