import os
import urllib.request
import logging
import argparse

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)

def download_webin_cli(version, dest_dir):
    jar_url = f"https://github.com/enasequence/webin-cli/releases/download/{version}/webin-cli-{version}.jar"
    dest_path = os.path.join(dest_dir, "webin-cli.jar")

    try:
        logging.info(f"Downloading webin-cli.jar version {version} to {dest_path} ...")
        urllib.request.urlretrieve(jar_url, dest_path)
    except Exception as e:
        logging.error(f"Failed to download webin-cli.jar: {e}")


def main():

    parser = argparse.ArgumentParser(description="Download the ENA webin-cli executable file")
    parser.add_argument("-d", "--dest-dir", default=".", help="Destination directory to save webin-cli.jar")
    parser.add_argument("-v", "--version", required=True, help="Release version/tag to download (default: latest)")
    args = parser.parse_args()
    download_webin_cli(args.version, args.dest_dir)


if __name__ == "__main__":
    main()
