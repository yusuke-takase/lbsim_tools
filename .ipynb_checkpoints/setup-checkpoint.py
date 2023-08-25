from setuptools import setup, find_packages

setup(
    name='lbsim_tools',  # パッケージ名（pip listで表示される）
    version="0.1.0",  # バージョン
    description="Packge for support the LiteBIRD simulation",  # 説明
    author='Yusuke Takase',  # 作者名
    packages=find_packages(),  # 使うモジュール一覧を指定する
    license='MIT'  # ライセンス
)