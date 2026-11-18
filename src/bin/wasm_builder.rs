//! Automated WebAssembly builder and GitHub Pages (`gh-pages`) deployer for `particle-accelarator`.

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};

const WASM_TARGET: &str = "wasm32-unknown-unknown";
const WASM_NAME: &str = "particle_accelarator.wasm";

fn run_cmd(program: &str, args: &[&str], cwd: Option<&Path>) -> Result<(), String> {
    let mut cmd = Command::new(program);
    cmd.args(args);
    if let Some(dir) = cwd {
        cmd.current_dir(dir);
    }
    let status = cmd
        .status()
        .map_err(|e| format!("Failed to spawn `{program}`: {e}"))?;
    if !status.success() {
        return Err(format!(
            "`{program} {}` exited with status {status}",
            args.join(" ")
        ));
    }
    Ok(())
}

fn capture_cmd(program: &str, args: &[&str], cwd: Option<&Path>) -> Result<String, String> {
    let mut cmd = Command::new(program);
    cmd.args(args);
    if let Some(dir) = cwd {
        cmd.current_dir(dir);
    }
    let output = cmd
        .output()
        .map_err(|e| format!("Failed to spawn `{program}`: {e}"))?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(format!(
            "`{program} {}` failed ({}): {stderr}",
            args.join(" "),
            output.status
        ));
    }
    Ok(String::from_utf8_lossy(&output.stdout).trim().to_string())
}

fn ensure_wasm_target() -> Result<(), String> {
    let installed = capture_cmd("rustup", &["target", "list", "--installed"], None)?;
    if !installed.lines().any(|line| line.trim() == WASM_TARGET) {
        println!("[wasm-builder] Installing Rust target `{WASM_TARGET}`...");
        run_cmd("rustup", &["target", "add", WASM_TARGET], None)?;
    } else {
        println!("[wasm-builder] Target `{WASM_TARGET}` already installed.");
    }
    Ok(())
}

fn validate_wasm_binary(wasm_path: &Path) -> Result<usize, String> {
    let bytes = fs::read(wasm_path)
        .map_err(|e| format!("Unable to read `{}`: {e}", wasm_path.display()))?;
    if bytes.len() < 8 || &bytes[0..8] != b"\x00asm\x01\x00\x00\x00" {
        return Err(format!(
            "Invalid WebAssembly magic header in `{}`",
            wasm_path.display()
        ));
    }
    Ok(bytes.len())
}

fn copy_file(src: &Path, dst: &Path) -> Result<(), String> {
    if let Some(parent) = dst.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("Failed to create dir `{}`: {e}", parent.display()))?;
    }
    fs::copy(src, dst)
        .map_err(|e| format!("Failed copying `{}` -> `{}`: {e}", src.display(), dst.display()))?;
    Ok(())
}

fn build_dist(root: &Path) -> Result<PathBuf, String> {
    ensure_wasm_target()?;

    println!("[wasm-builder] Compiling release WebAssembly library...");
    run_cmd(
        "cargo",
        &["build", "--target", WASM_TARGET, "--release", "--lib"],
        Some(root),
    )?;

    let wasm_src = root
        .join("target")
        .join(WASM_TARGET)
        .join("release")
        .join(WASM_NAME);
    let wasm_bytes = validate_wasm_binary(&wasm_src)?;
    println!(
        "[wasm-builder] Validated `{}` ({:.2} KiB).",
        wasm_src.display(),
        (wasm_bytes as f64) / 1024.0
    );

    let dist_dir = root.join("dist");
    if dist_dir.exists() {
        fs::remove_dir_all(&dist_dir)
            .map_err(|e| format!("Failed cleaning `dist/`: {e}"))?;
    }
    fs::create_dir_all(&dist_dir).map_err(|e| format!("Failed creating `dist/`: {e}"))?;

    copy_file(&wasm_src, &dist_dir.join(WASM_NAME))?;

    for asset in ["index.html", "styles.css", "app.js"] {
        copy_file(&root.join("web").join(asset), &dist_dir.join(asset))?;
    }

    for doc in [
        "SDLC.md",
        "SRS.md",
        "AGILE_TDD_PLAN.md",
        "SDD.md",
    ] {
        copy_file(
            &root.join("docs").join(doc),
            &dist_dir.join("docs").join(doc),
        )?;
    }

    let license_src = root.join("LICENSE");
    if license_src.exists() {
        copy_file(&license_src, &dist_dir.join("LICENSE"))?;
    }

    fs::write(dist_dir.join(".nojekyll"), b"")
        .map_err(|e| format!("Failed writing `.nojekyll`: {e}"))?;

    println!("[wasm-builder] Successfully assembled `dist/` bundle.");
    Ok(dist_dir)
}

fn deploy_gh_pages(root: &Path, dist_dir: &Path) -> Result<(), String> {
    let remote_url = capture_cmd("git", &["remote", "get-url", "origin"], Some(root))?;
    println!("[wasm-builder] Deploying `dist/` to `gh-pages` on `{remote_url}`...");

    let git_dir = dist_dir.join(".git");
    if git_dir.exists() {
        fs::remove_dir_all(&git_dir)
            .map_err(|e| format!("Failed removing old dist/.git: {e}"))?;
    }

    run_cmd("git", &["init"], Some(dist_dir))?;
    run_cmd("git", &["checkout", "-B", "gh-pages"], Some(dist_dir))?;
    run_cmd("git", &["add", "-A"], Some(dist_dir))?;
    run_cmd(
        "git",
        &[
            "commit",
            "-m",
            "deploy: publish Rust WASM particle accelerator simulator to gh-pages",
        ],
        Some(dist_dir),
    )?;
    run_cmd(
        "git",
        &["push", "-f", &remote_url, "gh-pages:gh-pages"],
        Some(dist_dir),
    )?;

    if git_dir.exists() {
        let _ = fs::remove_dir_all(&git_dir);
    }

    println!("[wasm-builder] Successfully deployed `gh-pages` branch!");
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = env::args().collect();
    let deploy = args.iter().any(|a| a == "--deploy");
    let root = env::current_dir().unwrap_or_else(|_| PathBuf::from("."));

    match build_dist(&root) {
        Ok(dist_dir) => {
            if deploy {
                if let Err(err) = deploy_gh_pages(&root, &dist_dir) {
                    eprintln!("[wasm-builder] Deployment error: {err}");
                    return ExitCode::FAILURE;
                }
            }
            ExitCode::SUCCESS
        }
        Err(err) => {
            eprintln!("[wasm-builder] Build error: {err}");
            ExitCode::FAILURE
        }
    }
}
