# GitHub'a yükleme rehberi

Bu arşiv GitHub'da yeni bir repository açmaya hazırdır; ZIP dosyasının kendisini yüklemek kod deposu oluşturmaz. Önce ZIP'i açıp içindeki `two-stage-distribution-network-control` klasörünü kullanın.

## 1. Başlangıç

GitHub'da New repository sayfasından `two-stage-distribution-network-control` adlı boş depo oluşturun. README, .gitignore veya LICENSE'i web arayüzünden otomatik eklemeyin: README ve .gitignore pakette bulunuyor. Ortak yazarların kod/veri dağıtım onayı ile üçüncü taraf dosyaların lisanslarını doğrulamadan depo herkese açık yapılmamalıdır. Önce private oluşturulabilir. Hakemlerin private depoyu otomatik açamayacağını unutmayın; bu durumda dergiye supplementary ZIP ayrıca yüklenmelidir.

## 2. Ubuntu/Linux üzerinde gönderim

```bash
unzip Alexis_IEEE_Access_GitHub_Repository_READY_20261022.zip
cd two-stage-distribution-network-control
python3 tools/verify_repository.py
git init
git add .
git commit -m "Add research code, fixed datasets and verification records"
git branch -M main
git remote add origin https://github.com/GITHUB_KULLANICI_ADI/two-stage-distribution-network-control.git
git push -u origin main
```

`GITHUB_KULLANICI_ADI` ifadesini kendi kullanıcı adınız veya kurum hesabınızla değiştirin. GitHub kimlik doğrulaması HTTPS için uygun token veya SSH anahtarı ile yapılır; şifrenizi kaynak koda, terminal geçmişine ya da bu sohbetin içine yapıştırmayın. `git commit` adı/e-posta isterse hesap ayarınızı yerel makinede yapılandırabilirsiniz.

## 3. Yayın sürümü

Son durumda bir sürüm etiketi (örneğin `v0.1.0`) tanımlanabilir. Makale yayımlanınca `CITATION.cff` içindeki makale bilgilerini ve DOI alanını gerçek yayın künyesine göre güncelleyin. Versiyonlu, kalıcı erişim gerekiyorsa GitHub sürümünü Zenodo'ya arşivlemek ayrı bir tercihtir.

## 4. Editör için bağlantı

GitHub bağlantısını Code and Data Availability açıklamasına ve uygun hakem yanıtına yalnızca depo gerçekten erişilebilir ve içerik doğrulanmış olduğunda ekleyin. URL henüz belli olmadığı için kaynak paketinde bağlantı uydurulmamıştır. Yayınlanmamış hakem yazışmalarını bu depoya yüklemeyin.
