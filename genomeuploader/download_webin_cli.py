import os
import sys
import urllib.request


def download_webin_cli(version, dest_dir):
    jar_url = f"https://github.com/enasequence/webin-cli/releases/download/{version}/webin-cli-{version}.jar"
    print(jar_url)
    dest_path = os.path.join(dest_dir, "webin-cli.jar")

    try:
        print(f"Downloading webin-cli.jar version {version} to {dest_path} ...")
        urllib.request.urlretrieve(jar_url, dest_path)
        print("Download complete")
    except Exception as e:
        print(f"Failed to download webin-cli.jar: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Download the ENA webin-cli JAR")
    parser.add_argument("-d", "--dest-dir", default=".", help="Destination directory to save webin-cli.jar")
    parser.add_argument("-v", "--version", required=True, help="Release version/tag to download (default: latest)")
    args = parser.parse_args()
    download_webin_cli(args.version, args.dest_dir)


if __name__ == "__main__":
    main()
